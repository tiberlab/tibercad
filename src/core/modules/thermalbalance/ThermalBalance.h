// $Id$

#ifndef _THERMALBALANCE_H_
#define _THERMALBALANCE_H_

#include "SimulationInterface.h"
#include "ElementSide.h"
#include "TiberLinearSystem.h"
#include "SimulationEnvironment.h"


/*!
 * 
 * \brief This is an example implementation of the MyPoisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class ThermalBalance : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~ThermalBalance(void);

    //! We need a public static creator function
    static ThermalBalance* create(const ModelOptions& options);



  protected:

    //! The initialization
    virtual void do_init(void);

    //! Parse the options from the input file
    virtual void parse_options(void);

    //! Setup the available variables
    virtual void do_setup_solution_variables(void);

    //! Solve the MyPoisson equation
    virtual void do_solve(void);

    //! Print some useful information
    virtual void do_print_info(void);

    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
					   const Material* mat) const;
   /*! \copydoc SimulationInterface::do_get_solution_vector() */
    virtual NumericVector<double>& do_get_solution_vector(void);


    //! We need to create boundary condition model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const Material* material_A, const Material* material_B) const;


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);


  //! Order the solution in correct mode
  virtual void build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results,
				       std::vector<std::string>& legend);


  private:


  typedef  std::map<const ElementSide, std::vector<double> >  SideData;

  


  SideData SD;
  std::vector<unsigned short int> node_conn;

  void compact(void);

  void do_partition_bis(void);

  void do_partition(void);

  bool is_on_GF_boundary(ElementSide elside);

  bool is_on_any_boundary(ElementSide elside);

  //Fourier Solution
  bool get_fourier_boundary(ElementSide elside,double a, double b, double c);

  double energy_conservation_check();

  double energy_conservation_check_traditional();

  bool is_fourier_solved;

  bool is_gray_solved;

  //------------------Gray solution--------------------------
  ID vec_spec;

  Point IntDir;

  Point dir;

  double d_omega;

  void from_nodal_to_cell(void);

  double get_boundary_value(ElementSide elside);
  ID gray_sys_number;
  std::vector< NumericVector<Number>* > sol_dir;
  std::vector< NumericVector<Number>* > thermal_flux;
  std::vector< NumericVector<Number>* > thermal_flux_nodal;
  NumericVector<Number>*  equilibrium_energy;
  NumericVector<Number>*  initial_energy;

    //!A Class that handle the angular integration
 class AngularIntegrator
 {
 public:

   // AngluarIntegrator();
   //~AngluarIntegrator(){};

   void compute_directions(void);

   void compute_custom_direction(std::vector<Point> custom_dir);

   void print_info(void);

   void print_info(ID k);
   
   ID dim;

   std::vector<double> d_omega;

   std::vector<double> theta_vec;

   std::vector<double> phi_vec;

   std::vector<Point> directions;

   std::vector<Point> dir;

   std::vector<Point> integrate_directions;

   std::vector<ID> spec;

   unsigned int theta_slices;

   unsigned int phi_slices;

   double d_theta;

   double d_phi;

   double total_angle;

   double weight;

   ID n_slices;
 };

 AngularIntegrator AngInt;

 struct options
  {
    ID custom;
    double max_iter;
    double max_error; //!< Max tollerance for self-consistent loop
    double ref_temp;
    double eq_temp;
    int DG;
    double work_units; //!< SI units, has to be consistent with the database parameters
    double initial_value;
    std::string macro_sim;
    std::vector<int> custom_dir;
   
    //New
    double equilibrium_energy;
    std::string first_guess;

    std::vector<Point> cd;

    bool diffusive;
    double s_0;
    double t_0;
    ID theta_slices;
    ID phi_slices;
   
    std::string partitioning;
    double threshold_value;
  };

 options myopts;

  //---------------------------------------------------------

  ID dim;

  void clear_system(const std::string& system_name);

  void from_cell_to_nodal(void);

  TiberLinearSystem*  system_fourier;

  void solve_fourier(void);

  void solve_gray(void);

  void do_init_fourier(void);

  void do_init_gray(void);

  struct thermal_options
  {
    
    bool automatic_partitioning;
    double threshold_value;
    double ms_error;
    ID ms_iter;
    bool fourier_guess;
    bool do_fourier;

  };
  thermal_options  opts;

  //  typedef  set<const Elem*>::iterator FourierIteratorFourier;

  //typedef  set<const Elem*>::const_iterator ConstIteratorFourier;

  std::set<const Elem*> FourierDomain;

  std::map<const Elem*, const Elem*> domain_boundary;

  std::set<ElementSide> BoundarySide;

  std::set<const Elem*> Domain;

  std::set<const Elem*> GlobalDomain;
  
  std::set<const Elem*> GrayDomain;

    //! These are the known solution variables
    /*!
     * This is an enum, but we use the string representation of 
     * the enum values to refer to solutions for plotting or 
     * for data exchange with other modules.
     *
     * \note Do \em not use (\c INVALID_ID - 1) or the strings \c RegionIDs
     * or \c materials as they are used to plot the materials/region IDs.
     *
     * \note The name "all" is used to plot all solutions
     */

    enum Solutions
    {
      temperature,       /*!< the Lattice Temperature */
      FourierTemp,       /*!< the Lattice Temperature */
      ThermalFlux,              /*!< the thermal flux */
      HeatSource,                /*!< the HeatSource */
      SolDir,
      Partition,
      ThermCond,
      thermal,
      DomainTest
    };

  //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    ThermalBalance(const ModelOptions& options);

    //! The assembly function
    static void assemble_fourier(EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_fourier(EquationSystems& es, const std::string& system_name);

     //! The assembly function
    static void assemble_gray(EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_gray(EquationSystems& es, const std::string& system_name);

     //! The assembly function
    static void assemble_global(EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_global(EquationSystems& es, const std::string& system_name);


    //! A static pointer to this
    static ThermalBalance* _this;


 

};

inline
bool
ThermalBalance::is_on_GF_boundary(ElementSide elside)
{

  //  bool is_on_GF_boundary = false;
  //if (domain_boundary.count(elside.elem()))
  // if (elside.elem()->neighbor(elside.side()) ==  domain_boundary[elside.elem()])
  //  is_on_GF_boundary = true;
  //----------------------------------	 
  //return is_on_GF_boundary;

  return BoundarySide.count(elside);

}

inline
bool
ThermalBalance::is_on_any_boundary(ElementSide elside)
{

  SimulationEnvironment& se = get_environment();
  //Check if Itnernal---------------------
  bool is_on_GF_boundary = false;
  if (domain_boundary.count(elside.elem()))
    if (elside.elem()->neighbor(elside.side()) ==  domain_boundary[elside.elem()])
      is_on_GF_boundary = true;
  //----------------------------------	   
  
  return se.is_outer_boundary(elside) || is_on_GF_boundary;

}

inline
void
ThermalBalance::assemble_fourier(EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_fourier(es, system_name);
}

inline
void
ThermalBalance::assemble_gray(EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_gray(es, system_name);
}

inline
void
ThermalBalance::assemble_global(EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_global(es, system_name);
}
#endif // _MYPOISSON_H_
