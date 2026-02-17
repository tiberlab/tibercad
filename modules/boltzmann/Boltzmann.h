/*  
 * This file is part of the tiberCAD module boltzmann.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file Boltzmann.h
 * \brief tiberCAD boltzmann module header.
 *
 * \note This file is part of module boltzmann.
 */


#ifndef TC_BOLTZMANN_H
#define TC_BOLTZMANN_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/geom/ElementSide.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/libMeshDefs.h"

class TiberLinearSystem;

/*!
 *
 * \brief This is an example implementation of the MyPoisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class TBDLLOCAL Boltzmann : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Boltzmann(void);

    //! We need a public static creator function
    static Boltzmann* create(const ModelOptions& options);



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
    virtual libMesh::NumericVector<double>& do_get_solution_vector(void);

    //! We need to create boundary condition model
    PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const;


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

    //! Get a mesh independent solution variable
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);

  //! Order the solution in correct mode
  virtual void build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results,
				       std::vector<std::string>& legend){};


  private:

  //! Check if the gray simulation should be solved
  bool is_gray;

   //! Variable for the boltzmann systems

  std::vector<ID> t_var;


  //! The index for the solid angle discretization
  ID solid_angle_iter;

  typedef  std::map<const ElementSide, std::vector<double> >  SideData;

  void get_gray_options(void);

  SideData SD;
  SideData SD_old;
  std::vector<unsigned short int> node_conn;

  //! Compute porosity
  double compute_porosity();

  //! Compute porosity
  double compute_view_factor(std::string S1, std::string S2);

  //! Compute surface to volume ratio
  double surface_to_volume_ratio();

  double compute_power_dissipated();

  double compute_effective_thermal_conductivity();

  double compute_effective_thermal_conductivity_elemental();

  double compute_power_emitted();

  void do_partition(void);

  bool is_on_GF_boundary(ElementSide elside);

  bool is_on_any_boundary(ElementSide elside);

  //Fourier Solution
  bool get_fourier_boundary(ElementSide elside,double a, double b, double c);

  double energy_conservation_check();

  double energy_conservation_check_traditional();

  bool is_fourier_solved;

  bool is_gray_solved;

  bool first_guess;

  //------------------Gray solution--------------------------
  ID vec_spec;

  Point IntDir;

  Point dir;

  double d_omega;

  void from_nodal_to_cell(void);

  double get_boundary_value(ElementSide elside, Point normal);
  ID gray_sys_number;
  std::vector< libMesh::NumericVector<Number>* > sol_dir;
  std::vector< libMesh::NumericVector<Number>* > thermal_flux;
  std::vector< libMesh::NumericVector<Number>* > thermal_flux_nodal;
  libMesh::NumericVector<Number>*  equilibrium_energy;

  libMesh::NumericVector<Number>*  initial_energy;

    //!A Class that handle the angular integration
 class AngularIntegrator
 {
 public:

   // AngluarIntegrator();
   //~AngluarIntegrator(){};

   void compute_directions(void);

   void compute_directions_bis(void);

   void compute_custom_direction(std::vector<Point> custom_dir);

   void print_info(void);

   void print_info(ID k);

   ID dim;

   double phi_zero;

   std::vector<double> d_omega;

   std::vector<double> theta_vec;

   std::vector<double> phi_vec;

   std::vector<Point> directions;

   std::vector<Point> dir;

   std::vector<ID> spec;

   unsigned int theta_slices;

   unsigned int phi_slices;

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
    double phi_zero;
    libMesh::RealGradient dist;
    //New
    double equilibrium_energy;
    std::string first_guess;

    bool compute_kappa;
    std::vector<Point> cd;
    std::string hot_contact;
    std::string cold_contact;
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


  void solve_boltzmann(void);

  void do_init_fourier(void);

  void do_init_gray(void);

  void do_init_boltzmann(void);

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

  //! Maximum temperature
  Real _Tmax;

  //! Minimum temperature
  Real _Tmin;


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
      LatticeTemp,       /*!< the Lattice Temperature */
      NormalizedLatticeTemp,       /*!< the Normalized Lattice Temperature */
      NormalizedThermalFlux,       /*!< the Normalized ThermalFlux */
      FourierTemp,       /*!< the Lattice Temperature */
      ThermalFlux,              /*!< the thermal flux */
      HeatSource,                /*!< the HeatSource */
      SolDir,
      Partition,
      ThermCond,
      EffectiveKappa,
      thermal,
      DomainTest,
      GRAY,
      MaxTemp
    };

  //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    Boltzmann(const ModelOptions& options);

    //! The assembly function
    static void assemble_fourier(libMesh::EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_fourier(libMesh::EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_boltzmann(libMesh::EquationSystems& es, const std::string& system_name);

    //! The assembly function
    static void assemble_boltzmann(libMesh::EquationSystems& es, const std::string& system_name);

     //! The assembly function
    static void assemble_gray(libMesh::EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_gray(libMesh::EquationSystems& es, const std::string& system_name);

     //! The assembly function
    static void assemble_global(libMesh::EquationSystems& es, const std::string& system_name);

    //! The real assembly function
    void do_assemble_global(libMesh::EquationSystems& es, const std::string& system_name);


    //! A static pointer to this
    static Boltzmann* _this;




};

inline
bool
Boltzmann::is_on_GF_boundary(ElementSide elside)
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
Boltzmann::is_on_any_boundary(ElementSide elside)
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
Boltzmann::assemble_fourier(libMesh::EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_fourier(es, system_name);
}

inline
void
Boltzmann::assemble_gray(libMesh::EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_gray(es, system_name);
}


inline
void
Boltzmann::assemble_boltzmann(libMesh::EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_boltzmann(es, system_name);
}

inline
void
Boltzmann::assemble_global(libMesh::EquationSystems& es, const std::string& system_name)
{
  _this->do_assemble_global(es, system_name);
}
#endif // TC_MYPOISSON_H
