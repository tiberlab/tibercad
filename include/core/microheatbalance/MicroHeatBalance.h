// $Id$

#ifndef _MICROHEATBALANCE_H_
#define _MICROHEATBALANCE_H_



//------------------------------------------------------------------------------


#include "SimulationInterface.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"
#include "HeatModel.h"

#include "elem.h"


class TiberLinearSystem;
class Device;
class Mesh;


//!  Class to solve heat transport problem
class MicroHeatBalance : public SimulationInterface
{

 public:


   enum Variables
   {
     UNKNOWN = 0,
     TEMPERATURE,
     JQX,
     JQY,
     JQZ,
     E0
   };

  //!options that we need for this simulation
  struct options
  {
    ID alternative;
    double max_iter;
    double max_error; //!< Max tollerance for self-consistent loop
    double ref_temp;
    double eq_temp;
    double work_units; //!< SI units, has to be consistent with the database parameters
    double initial_value;
    std::string macro_sim;
     libMeshEnums::QuadratureType  quadrature_type;
    /**
     * The order of gauss integration
     */
    libMeshEnums::Order integration_order;
    double scale;

    //New
    double equilibrium_energy;
    std::string first_guess;
  
  };


  /*!
   * \copydoc SimulationInterface::get_solution_secure(const Elem*,
   * const std::vector<Point>&, const std::vector<ID>&,
   * std::vector<std::vector<double> >&)
   */
  virtual void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);


  /*!
   * \copydoc SimulationInterface::get_solution_secure(const Elem*,
   * const std::vector<ID>&, std::vector<std::vector<double> >&)
   */
  virtual void get_solution_secure(const Elem* elem,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);

  virtual  NumericVector< double > & do_get_solution_vector(void);

  //!Constructor
  MicroHeatBalance(const ModelOptions& options);

  //!Destructor
  virtual ~MicroHeatBalance();

  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);


  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const
    throw (ModelErrorException);

  static void assemble_heat_matrix(EquationSystems& es,
				   const std::string& system_name);

  static void assemble_macro_heat_matrix(EquationSystems& es,
				  const std::string& system_name);

  //!Create an MacroHeatBalance object
  static MicroHeatBalance*  create(const ModelOptions& options);

  /*! \copydoc SimulationInterface::build_integrated_quantities() */
  virtual void build_integrated_quantities(std::vector<double>& values);


  /*! \copydoc SimulationInterface::build_integrated_quantities_description()
   */
  virtual void build_integrated_quantities_description(
      std::vector<std::string>& legend,
      std::vector<std::string>& description);

   /*! \copydoc SimulationInterface::do_print_info() */
    virtual void do_print_info(void);

 private:

  double s_0;
  double t_0; 
  //Boundary values
  std::map<const Node*, std::vector<double> > bv;

  typedef  std::map<const Node*, std::vector<double> >  BoundaryData;

  void compute_fourier_solution(void);

  double reference_temperature;

  std::map< const Elem*,double >  eq_energy;

  //std::map< const Elem*,RealGradient> thermal_flux;

  //  double scale;

  double vg;
  double tg;
  double cg;

  std::string heat_legend;

  std::set<ID> JQ_var;

  std::set<ID> e0_var;

  void compute_flux(void);

  void reset_flux(void);

   //! Quadrature point along the face of the element
  const std::vector<Point> qface_point;

  //! A pointer to heat simulation
  MicroHeatBalance* _heat_simul;

  //EquationSystems *equation_systems;

  HeatModel* heat_model;

  void init_heat_model(const Elem* elem);


  std::string system_name;

  TiberLinearSystem* my_system;

  TiberLinearSystem*  fourier_system;
  //! Order the solution in correct mode
  virtual void 	build_nodal_results(const std::set< std::string > &variables,
				     std::vector< double > &results,
				     std::vector< std::string > &legend);


  //! Order the solution in correct mode
  virtual void build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results,
				       std::vector<std::string>& legend){};

  //! Calculate Power Dissipated
  /*!
   * Integrates numerically over the boundary elements.
   * The power dissipated is then:
   *
   * \f[P = \int_{\Gamma} -\kappa \nabla T \cdot \mathbf{N} \mathrm{d}\Gamma \f]
   */

  // void compute_total_flux(void);

  //void compute_total_equilibrium_energy(void);

  double calculate_power_emitted(void);

  // double calculate_power_dissipated_rstf(void);

  void calculate_power_dissipated(double& power_dissipated, double& error);

  static MicroHeatBalance* static_this;

  options myopts;

  //!non-static method that actually does matrix assembling
  void do_assemble(EquationSystems& es, const std::string& system_name);

  void do_macro_assemble(EquationSystems& es, const std::string& system_name);

  static Device* _device;

  //!Dimension of meshmap
  short dim;
  
  ID vec_spec;

  Point IntDir;

  Point dir;

  double d_omega;

  //!Pointer to mesh
  //Mesh* mesh;
  const MeshBase* mesh;
  //! Object that handles angular integration

  //std::map<ID,AutoPtr<NumericVector<Number> > > dir_solution;

 //  AutoPtr<NumericVector<Number> > thermal_flux_x;
//   AutoPtr<NumericVector<Number> > thermal_flux_y;
//   AutoPtr<NumericVector<Number> > thermal_flux_z;


  

  ID dir_iter;

  //!A Class that handle the angular integration
 class AngularIntegrator
 {
 public:

   // AngluarIntegrator();
   //~AngluarIntegrator(){};

   void compute_directions(void);

   void compute_alternative_directions(void);

   void compute_very_alternative_directions(void);

   void compute_very_alternative_directions2(void);

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

 void compute_equilibrium_energy(void);

  NumericVector<Number>* equilibrium_energy;

  NumericVector<Number>* equilibrium_energy_new;

  //std::vector< AutoPtr<NumericVector<Number> > > sol_dir(4);


  AngularIntegrator AngInt;

  std::vector< NumericVector<Number>* > sol_dir;
  std::vector< NumericVector<Number>* > sol_dir_new;
  std::vector< NumericVector<Number>* > thermal_flux;

 protected:

    /*! \copydoc SimulationInterface::convert_variable_name_to_id() */
  virtual ID convert_variable_name_to_id(const std::string& variable_name) const;

  //! \copydoc  SimulationInterface::do_init()
  virtual void 	do_init (void);

  //!Do the solve
  virtual void do_solve (void);

  //!Parse the options
  virtual void 	parse_options(void);


};

#endif
