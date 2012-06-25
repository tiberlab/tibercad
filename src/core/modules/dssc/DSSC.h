// $Id$

#ifndef _DSSC_H_
#define _DSSC_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "Device.h"
#include "SolverException.h"

// Libmesh includes
#include "libmesh_common.h"
#include "enum_order.h"


// C++ includes
#include <vector>
#include <set>
#include <map>


// forward declarations
//class Device;
class Boundary;
class MeshBase;
class Elem;
class Point;
class Node;
class EquationSystems;


template<typename T> class DenseMatrix;
template<typename T> class NumericVector;
template<typename T> class SparseMatrix;


//! The main class to perform DSSC calculations
class DSSC : public SimulationInterface
{

  public:

  public:

    //! The variables that can be provided
    enum Solutions
    {
      //Ec,             /*!< conduction band edge */
      //Ev,             /*!< valence band edge */
      eQFermi,          /*!< electron electro-chemical potential */
      IQFermi,          /*!< iodide electro-chemical potential */
      I3QFermi,         /*!< triiodide electro-chemical potential */
      CQFermi,           /*!< cation electro-chemical potential */
      ElPotential,      /*!< electric potential */
      Eredox,           /*!< Eredox potential > */
      //Ec0,              /*!< bare conduction band edge */
      //Ev0,              /*!< bare valence band edge */
      //Eg,               /*!< band gap */
      //ConductionBands,  /*!< the conduction band energies */
      //ValenceBands,     /*!< the valence band energies */
      eDensity,         /*!< electron density */
      IDensity,         /*!< iodide density */
      I3Density,         /*!< triiodide density */
      CDensity,         /*!< cation density */
      Traps,            /*!< trap density */
      eMobility,        /*!< electron mobility */
      //eConductivity,    /*!< electron conductivity */
      ElField,          /*!< electric field vector */
      CurrentDensity,   /*!< total electric current density */
      //eFlux,            /*!< electron flux */
      //hFlux,            /*!< hole flux */
      eCurrentDensity,  /*!< electron current density */
      ICurrentDensity,  /*!< iodide current density */
      I3CurrentDensity,  /*!< triiodide current density */
      CCurrentDensity,  /*!< cation current density */
      //IonizedElectronTraps, /*!< trapped electron density */
      Generation,        /*!< generation term > */
      NetRecombination = 100,  /*!< base number for recombination models */
      ContactCurrent   = 200,  /*!< base number for contact currents */
      ContactVoltage   = 300   /*!< base number for contact voltages */
    };



    //! The variables that can be provided
    //enum Variables
    //{
    //  UNKNOWN = 0,
    //  ELPOTENTIAL,      /*!< electric potential */
    //  QFERMIE           /*!< electron electro-chemical potential */
    //};


    //! Constructo(const ModelOptions& options);
    DSSC(const ModelOptions& options);

    //! Destructor
    virtual ~DSSC(void);


    //! Create an DSSC object
    static DSSC* create(const ModelOptions& options);


    /*! \copydoc SimulationInterface::create_physical_model() */
    virtual PhysicalModel*
      create_bulk_model(const ModelOptions& options,
          const Material* mat) const;


    /*! \copydoc SimulationInterface::create_boundary_model() */
    virtual BoundaryProperties*
      create_boundary_model(const ModelOptions& options) const
      throw (ModelErrorException);


    //! Get the mesh
    /*!
     * \return a constant reference to the simulation mesh
     */
    MeshBase& get_mesh(void) const;


    //! Makes a first guess of the equilibrium potential
    /**
     * It sets every node to its equilibrium potential.
     */
    void guess_equilibrium(void);


    //! Get the number of nonlinear iterations needed for the solution
    unsigned int get_n_nonlinear_iterations(void) const;


    //! Get the final residual norm of the solution
    double get_final_residual(void) const;


    //! Get the boundary currents indexed by boundary descriptor
    const std::map<const Boundary*, double>&
      get_boundary_currents(void) const;




  protected:


    //! Initialize the equation system
    /*!
     * This has to be called before any call to methods which access the
     * equation system object. So do it early.
     */
    virtual void do_init(void);


    /*! \copydoc SimulationInterface::do_equilibrium() */
    virtual void do_equilibrium(void);


    //! Solve the drift-diffusion problem.
    /*!
     * If adaptive mesh refinement was enabled for this solver
     * run, it will be deactivated afterwards and has to be re-enabled
     * explicitly.
     */
    virtual void do_solve(void);


    /*! \copydoc SimulationInterface::do_print_info() */
    virtual void do_print_info(void);


    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);


    /*! \copydoc SimulationInterface::solution_vector() */
    virtual NumericVector<double>& solution_vector(void);


//    /*! \copydoc SimulationInterface::build_nodal_results() */
//    virtual void build_nodal_results(const std::set<std::string>& variables,
//        std::vector<double>& results, std::vector<std::string>& legend);


//    /*! \copydoc SimulationInterface::build_elemental_results() */
//    virtual void build_elemental_results(const std::set<std::string>& variables,
//        std::vector<double>& results, std::vector<std::string>& legend);

    //! Setup the available variables
    virtual void do_setup_solution_variables(void); 



    /*!
     * \copydoc SimulationInterface::get_solution_secure(const Elem*,
     *  std::map<ID, std::vector<double> >&, const std::vector<Point>&)
     */
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);


    //! We override this to not produce too many files ...
    virtual void plot_globaldata(void) {};


    /*!
     * \copydoc SimulationInterface::get_solution_secure(
     *  std::map<ID, std::vector<double> >&)
     */
    virtual void get_solution_secure(
        std::map<ID, std::vector<double> >& values);


    /*! \copydoc SimulationInterface::build_integrated_quantities() */
    //virtual void build_integrated_quantities(std::vector<double>& values);


    /*! \copydoc SimulationInterface::build_integrated_quantities_description()
     */
    //virtual void build_integrated_quantities_description(
    //    std::vector<std::string>& legend,
    //    std::vector<std::string>& description);


    /*! \copydoc SimulationInterface::do_get_solution_vector() */
    virtual NumericVector<double>& do_get_solution_vector(void);


    /*! \copydoc SimulationInterface::do_maximum_norm_of_difference() */
    virtual double do_maximum_norm_of_difference(ID id);


    /*! \copydoc SimulationInterface::convert_variable_name_to_id() */
    //virtual ID convert_variable_name_to_id(const std::string& variable_name) const;


    /*!
     * \copydoc SimulationInterface::get_solution_secure(const Elem*,
     * const std::vector<ID>&, std::vector<std::vector<double> >&)
     */
    //virtual void get_solution_secure(const Elem* elem,
    //    const std::set<ID>& ids,
    //    std::vector<std::map<ID, double> >& values);


    /*!
     * \copydoc SimulationInterface::get_solution_secure(const Elem*,
     * const std::vector<Point>&, const std::vector<ID>&,
     * std::vector<std::vector<double> >&)
     */
    //virtual void get_solution_secure(const Elem* elem,
    //    const std::vector<Point>& p,
    //    const std::set<ID>& ids,
    //    std::vector<std::map<ID, double> >& solution);


  private:


    // for nicer code
    typedef std::map<const Boundary*, double> ContactData;
    typedef std::map<const Node*, Boundary*> BoundaryNodeList;


    struct ConductivityScaling
    {
      double n;
      double I;
      double I3;
      double C;
    };


    ConductivityScaling _cond_scaling;


    //! A static reference to \c this
    /*!
     * This is needed during matrix assembly, which is a static method.
     */
    static DSSC* _this;


    //! An internal pointer to the device
    Device* _device;


    /*!
     * A list of nodes with dirichlet boundary conditions
     */
    BoundaryNodeList _dirichlet_nodes;


    //! The total number of cations
    double _cation_amount;


    //! The total amount of iodine
    double _iodine_amount;


    /*!
     * If @c true, the equation system needs to be rebuilt
     */
    bool _rebuild_eq_system;


    //! Do only Poisson if true
    bool _poisson_only;


    //! Whether we do EIS or not
    bool _do_EIS;

    //! The current frequency (it's here because it's a hack...)
    double _frequency;

    //! Do the actual EIS calculation
    void do_EIS(void);


    //! Tells if we are doing only Poisson
    bool poisson_only(void) const;


    //! coordinates of the contact from where the light comes 
    Point _x0;


    //! The type of scaling
    Scaling::ScalingType _scaling_type;


    //! The boundary currents
    /*!
     * Currents are calculated after each solve step.
     */
    ContactData _boundary_currents;
    
    

    //! The voltages of the former solve step
    ContactData _voltages;


    //! The local density scaling
    std::map<const Node*, std::vector<double> > local_scaling_;


    //! If true, local density scaling should be applied
    bool do_local_scaling_;


    /*!
     * \brief A list of nodes which lie on an inner boundary between
     * TiO2 and electrolyte
     */
    std::set<const Node*> _internal_boundary_nodes;


    //! Calculate the local density scaling on each node
    void build_local_scaling(void);


    /**
     * The number of nonlinear iterations needed
     */
    unsigned int _n_nonlinear_iterations;


    /**
     * The final residual norm
     */
    double _final_residual;


    //! disable the copy constructor
    DSSC(const DSSC& rhs);


    //! disable the copy assignment operator
    DSSC& operator=(const DSSC& rhs);


    //! Parse the options which will not change between calls to solve()
    void parse_const_options(void);


    //! Rebuild the equation system if needed
    void rebuild_equation_system(void);


    /*!
     * \brief Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(Scaling::ScalingType type = Scaling::UNITS);
    void compute_scaling_only(Scaling::ScalingType type = Scaling::UNITS);


    //! Fills the dirichlet nodes data structure.
    void find_dirichlet_nodes(void);


    //! Find nods on boundary TiO2/electrolyte
    void find_internal_boundary_nodes(void);


    //! Tells if node lies on an inner TiO2/electrolyte boundary
    bool is_internal_boundary_node(const Node* node) const;


    //! Reset solver environment.
    /*!
     * Deletes only the \p EquationSystems object, without
     * touching simulation voltages and solutions.
     */
    void reset_solver(void);


    //! Cleanup solver environment.
    /*!
     * Deletes the \p EquationSystems object, the simulation voltages
     * vector and the solution and variable vectors.
     */
    void cleanup_solver(void);


    //! Solve using an iterative Gummel scheme
    //void solve_gummel(void);


    //! Do a Newton type iteration
    void do_newton(void);


    //! Calculate the terminal currents
    /*!
     * Calls \c calculate_currents_surfint() or
     * calculate_currents_rstf()
     */
    void calculate_currents(void);


    //! Calculate terminal currents
    /*!
     * Uses the Ramo-Shockley test functions and integrates over the
     * volume.
     *
     * Assuming electron and hole generation-recombination terms to be
     * equal, one can write:
     * \f[\left(-\nabla(\mathbf{j}_n + \mathbf{j}_p), hl\right) =
     * -\int_\Omega h_l \nabla (\mathbf{j}_n + \mathbf{j}_p)\mathrm{d}V = 0\f]
     * where \f$h_l|_{\Gamma_j} = \delta_{lj}, h_l \in H^1\f$ is the test
     * function for the contact \it l
     *
     * Using Gauss in the scalar product, one gets
     * \f[ 0 = - \int_{\partial\Omega}h_l(\mathbf{j}_n +
     * \mathbf{j}_p)\mathrm{\mathbf{S}}
     * + \int_\Omega \nabla h_l (\mathbf{j}_n +
     * \mathbf{j}_p)\mathrm{d}V\f]
     * The first term of the left hand side is exactly the terminal current
     * of contact \it l due to our choice of the testfunction, therefore
     * \f[ I_l = \int_\Omega \nabla h_l (\mathbf{j}_n +
     * \mathbf{j}_p)\mathrm{d}V\f]
     *
     * In this implementation we choose the test function to be (where \it i
     * is a node index and \f$\psi_i\f$ the FEM basis function associated with
     * node \it i):
     * \f[h_l = \sum_{i \in \Gamma_l}\psi_i\f]
     */
    void calculate_currents_rstf(void);


    //! similar routine to compute complex current for Impedance spectroscopy
    void calculate_currents_rstf_EIS(std::map<const Boundary*, double>& curr_R,
        std::map<const Boundary*, double>& curr_I);


    //! Assemble the residual vector or the jacobian matrix
    /*!
     * This method gets called from the underlying nonlinear solver
     * library. It calls the real assembly routines.
     */
    static void assemble_system(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);


    //! Assembles the residual vector or the jacobian matrix
    /*!
     * Assembles the residual vector or the jacobian matrix for
     * the equation system with @c Coupling T.
     *
     * This implementation uses standard FEM.
     */
    void do_assembly(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);


    //! The assembly function for EIS
    static void assemble_EIS(EquationSystems& es, const std::string& system_name);


    //! Assmebly of the Impedance Spectroscopy Jacobian matrix
    void do_assembly_frequency(EquationSystems& es, const std::string& system_name);


    //! Find the open circuit potential and densities
    void get_OC_values(void);

};


//
// inline member functions
//


inline
DSSC*
DSSC::create(const ModelOptions& options)
{
  return new DSSC(options);
}


inline
unsigned int
DSSC::get_n_nonlinear_iterations(void) const
{
  return _n_nonlinear_iterations;
}


inline
double
DSSC::get_final_residual(void) const
{
  return _final_residual;
}


inline
const std::map<const Boundary*, double>&
DSSC::get_boundary_currents() const
{
  return _boundary_currents;
}



inline
MeshBase&
DSSC::get_mesh(void) const
{
  return _device->get_mesh();
}


inline
bool
DSSC::poisson_only(void) const
{
  return _poisson_only;
}


inline
bool
DSSC::is_internal_boundary_node(const Node* node) const
{
  bool result = false;
  if (_internal_boundary_nodes.find(node) != _internal_boundary_nodes.end())
    result = true;

  return result;
}


#endif // _DSSC_H_
