// $Id$

#ifndef _DRIFTDIFFUSION_H_
#define _DRIFTDIFFUSION_H_

#include "SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "Scaling.h"
#include "DDevice.h"

// Libmesh includes
#include "libmesh_common.h"
#include "enum_order.h"

// PETSc include
#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscksp.h>
}
#else
# include <petscksp.h>
#endif

// C++ includes
#include <vector>
#include <map>

// forward declarations
class DD::Device;
class BoundaryDescriptor;
class Mesh;
class Elem;
class Node;
class EquationSystems;
class NonlinearImplicitSystem;

template<typename T> class NumericVector;
template<typename T> class SparseMatrix;

template<typename T> class TiberPetscNonlinearSolver;
template<typename T> class NonlinearSolver;

//! The main class to perform standard drift-diffusion calculations
/*!
 * TODO
 * Some more details
 */
class DriftDiffusion
{
  public:
 
    //! The solver methods that can be used
    enum SolverMethod
    {
      NEWTON,
      GUMMEL
    };
    
    //! This class defines parameters used by the underlying
    //! nonlinear solver
    /*!
     *  For details refer to the
     * <A HREF="http://www-unix.mcs.anl.gov/petsc/petsc-2">online
     * PETSc documentation</A>.
     */
    class SolverParameters
    {
      public:
        SolverParameters(void);
        
        SolverParameters(const SolverParameters& rhs);
        
        SolverParameters& operator=(const SolverParameters& rhs);
      
        double nonlinear_tolerance;
        double nonlinear_abs_tolerance;
        unsigned int nonlinear_max_iterations;
        double linear_tolerance;
        double linear_abs_tolerance;
        unsigned int linear_max_iterations;
        
        //! The line search maximum step size per grid point
        double ls_maxstep;

        //! The linear (KSP) solver type
        /*!
         * This defines what type of linear solver to be used. For
         * details refer to the
         * <A HREF="http://www-unix.mcs.anl.gov/petsc/petsc-2">online
         * PETSc documentation</A>.
         *
         * \note
         * Usually \c KSPBCGS or \c KSPBCGSL seem to be the most stable
         * solver types
         */
        KSPType ksp_type;

        //! The preconditioner (PC)
        /*!
         * This defines the type of preconditioner to be used. For
         * details refer to the
         * <A HREF="http://www-unix.mcs.anl.gov/petsc/petsc-2">online
         * PETSc documentation</A>.
         *
         * \note
         * In most cases \c PCILU seems to be a good choice, but sometimes
         * there are problems with zero pivot values. In this case, the use of
         * \c PCJACOBI can solve it.
         */
        PCType pc_type;

      private:

        friend class DriftDiffusion;
    };
    
    /**
     * This class defines various parameters that control a
     * drift-diffusion calculation
     */
    class Options
    {
      public:

        Options(void);
        
        Options(const Options& rhs);

        Options& operator=(const Options& rhs);

        /**
         * Enable or disable mesh refinement
         */
        bool mesh_refinement;

        /**
         * Maximum number of refinement steps
         */
        int max_refinement_steps;

        /**
         * Maximum allowable refinement level
         *
         * The refinement hierarchy will not get deeper than this
         * value, i.e. an element will not be refined more than
         * @c max_refinement_level times.
         */
        int max_refinement_level;

        /**
         * The error fraction for refinement
         *
         * All elements with percentual error bigger than @c refine_fraction
         * will get refined.
         */
        float refine_fraction;

        /**
         * The error fraction for coarsening
         *
         * All elements with percentual error lower than @c coarsen_fraction
         * will get coarsened.
         */
        float coarsen_fraction;

        /**
         * The tolerance at which mesh refinement will stop
         */
        double refinement_tolerance;

        /**
         * The minimum allowable voltage step size
         */
        double min_voltage_step;

        /**
         * The order of gauss integration
         */
        libMeshEnums::Order integration_order;

        /**
         * The approximation order for the finite elements
         */
        libMeshEnums::Order approximation_order;

        /**
         * The solver method
         *
         * Can be @c NEWTON or @c GUMMEL
         */
        SolverMethod solver_method;

        /**
         * The maximum number of iteration steps for
         * the Gummel method
         */
        int max_gummel_iterations;

        /**
         * The nonlinear/linear (PETSc-)solver parameters
         */
        SolverParameters solver_params;

        /**
         * The units in which the mesh object is given
         */
        double mesh_units;

        /**
         * The type of scaling to be applied
         *
         * Can be one of @c NONE, @c UNITS or @c DEMARI
         */
        Scaling::ScalingType scaling_type;

        //! The type of coupling to be used
        /*!
         * We safe it as int, because we want to assign values by
         * using logic operators, e.g. \code ELECTRONS | POISSON \endcode
         */
        int coupling;

      private:
        
        /**
         * The maximum electron density
         */
        double n_max;
        
        /**
         * The maximum hole density
         */
        double p_max;

        /**
         * The density scaling factor for the electron current equation
         */
        double C0_e;

        /**
         * The density scaling factor for the hole current equation
         */
        double C0_h;

        friend class DriftDiffusion;
    };

      
    DriftDiffusion(DD::Device* device);

    DriftDiffusion(DD::Device* device, DriftDiffusion::Options& params);

    ~DriftDiffusion(void);
    
    /**
     * @returns a reference to the device to be solved
     */
    const DD::Device& get_device(void) const;
    
    /**
     * Set a new device to be simulated
     *
     * After a call to this function, the equation system object
     * and the simulation voltages map are cleared and have to be
     * rebuilt in the \p solve method.
     */
    void set_device(DD::Device* device);

    /**
     * @returns a reference to the simulation options
     */
    Options& get_options(void);

    /**
     * Enables adaptive mesh refinement.
     */
    void enable_mesh_refinement(void);

    /**
     * Disables adaptive mesh refinement.
     */
    void disable_mesh_refinement(void);

    /**
     * Set new simulation options.
     */
    void set_options(const Options& options);

    //! Get the mesh
    /*!
     * \return a constant reference to the simulation mesh
     */
    Mesh& get_mesh(void) const;
    
    /**
     * @returns a constant reference to the current scaling parameters
     */
    const Scaling& get_scaling(void) const;

    /**
     * Set the simulation voltage for the boundary named \p boundary.
     */
    // TODO throw exception if boundary non-existing
    void set_simulation_voltage(const std::string& boundary,
        double voltage);

    /**
     * Remember the current solution for future restart
     */
    void remember_current_solution(void);

    /**
     * Reset to the remembered solution
     */
    void set_to_remembered_solution(void);


    /**
     * Solve the drift-diffusion problem.
     *
     * If adaptive mesh refinement was enabled for this solver
     * run, it will be deactivated afterwards and has to be re-enabled
     * explicitly.
     */
    void solve(void);

    /**
     * Solve the drift-diffusion problem starting from the equilibrium
     * solution.
     *
     * If adaptive mesh refinement was enabled for this solver
     * run, it will be deactivated afterwards and has to be re-enabled
     * explicitly.
     */
    void solve(bool restart);
   
//    void solve(SolverMethod method);

    /**
     * @returns the variable names in a vector.
     */
    const std::vector<std::string>& get_variable_names(void) const;

    /**
     * @returns the current nodal solution vector.
     *
     * The vector is ordered as
     *   [ var1\@node1 ... varn\@node1 var1\@node2 ... varn\@node2 ... ]
     * where the ordering \p var1 ... \p var2 is the same as in the
     * method \p get_variable_names()
     */
    const std::vector<Number>& get_solution(void) const;

    /**
     * @returns the number of nonlinear iterations needed for the solution
     */
    unsigned int get_n_nonlinear_iterations(void) const;

    /**
     * @returns the final residual norm of the solution
     */
    double get_final_residual(void) const;

    /**
     * @returns the boundary currents indexed by boundary descriptor
     * pointers.
     */
    const std::map<const BoundaryDescriptor*, double>&
      get_boundary_currents(void) const;

    /**
     * @returns the integrated current at the artificial (non-contact)
     * boundary of the device, which should be zero.
     */
    double get_artificial_boundary_current(void);


    /**
     * fill @c densities with electron and hole densities
     *
     * TODO add other values
     */
    void build_densities(std::vector<double>& densities,
        std::vector<std::string>& names);


  private:

    /**
     * Are we doing equilibrium or nonequilibrium calculation?
     */
    enum CalculationType
    {
      EQUILIBRIUM,
      NONEQUILIBRIUM
    };

    // for nicer code
    typedef std::map<const BoundaryDescriptor*, double> ContactData;
    typedef std::map<const Node*, const BoundaryDescriptor*> BoundaryNodeList;
    typedef TiberPetscNonlinearSolver<Real> SolverClass;

    //! A static reference to \c this
    /*!
     * This is needed during matrix assembly, which is a static method.
     */
    static DriftDiffusion* _this;

    /**
     * The device to be solved by this DriftDiffusion object
     */
    DD::Device* _device;

    /**
     * A list of nodes with dirichlet boundary conditions
     */
    BoundaryNodeList _dirichlet_nodes;

    /**
     * The equation system for this device
     */
    EquationSystems* _eq_system;
    
    /**
     * If @c true, the equation system needs to be rebuilt
     */
    bool _rebuild_eq_system;

    /**
     * The @c Options to be used
     */
    Options _options;
    
    /**
     * The scaling parameters
     */
    Scaling _scaling;

    /**
     * The simulation voltages of the previous solve step
     */
    ContactData _old_sim_voltages;
    
    /**
     * The simulation voltages for the next simulation
     */
    ContactData _simulation_voltages;
    
    /**
     * The remembered simulation voltages
     */
    ContactData _remembered_voltages;

    /**
     * The boundary currents
     *
     * Currents are calculated after each solve step.
     */
    ContactData _boundary_currents;

    /**
     * The variable names in the same order as they appear in the
     * solution vector
     */
    std::vector<std::string> _variables;

    /**
     * The nodal solution in the order
     * [ var1\@node1 ... varn\@node1 var1\@node2 ... varn\@node2 ... ]
     */
    std::vector<Number> _solution;

    /**
     * The number of nonlinear iterations needed
     */
    unsigned int _n_nonlinear_iterations;

    /**
     * The final residual norm
     */
    double _final_residual;

    // disable the copy constructor and assignment operator
    DriftDiffusion(const DriftDiffusion& rhs);
    DriftDiffusion& operator=(const DriftDiffusion& rhs);

    /**
     * Set the options for the PETSc solver as given in @c SolverParameters
     */
    void set_solver_params(NonlinearSolver<Number>& solver,
        CalculationType calc_type);

    /**
     * Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(Scaling::ScalingType type = Scaling::UNITS);

    /**
     * prepare data structures used for solving a system
     */
    void prepare_solver(void);

    /**
     * Fills the dirichlet nodes data structure.
     */
    void find_dirichlet_nodes(void);

    /**
     * Sets the boundary values for each node on a Dirichlet type
     * boundary
     *
     * This function could be called after
     * \p calculate_next_sim_voltage() and after the static parameters
     * were set up.
     */
    //void set_dirichlet_values(void);

    /**
     * Reset solver environment.
     *
     * Deletes only the \p EquationSystems object, without
     * touching simulation voltages and solutions.
     */
    void reset_solver(void);

    /**
     * Cleanup solver environment.
     *
     * Deletes the \p EquationSystems object, the simulation voltages
     * vector and the solution and variable vectors.
     */
    void cleanup_solver(void);

    /**
     * Initializes the equation system \p system and prepares it
     * to be solved
     */
    void initialize_eq_system(EquationSystems& system);

    //! Get the equation system
    EquationSystems& get_equation_system(void);

    void solve_newton(bool restart);

    void solve_gummel(bool restart);

    void calculate_currents(void);


    /**
     * Calculates the next simulation point and returns
     * the voltage step
     */
    double calculate_new_simulation_voltages(void);


    /**
     * @returns the side number of the top level element, if the side
     * \p side of element \p elem is on a boundary, -1 if not.
     *
     * It is supposed that \p top_parent is the top parent of \p elem
     * and \p side lies on a device boundary
     */
    static int find_boundary(const Elem* elem, int side,
      const Elem* top_parent);

    //! Assign boundary value coefficients
    /*!
     * Assigns boundary value coefficient in a form
     *
     *   du/dn = - coeff * u + value
     *
     * NOTE: this method assumes mixed or von Neumann boundary conditions.
     */
    static void assign_boundary_values(double& coeff, double& value,
        const std::vector<double>& coefficients);

    //! Makes a first guess of the equilibrium potential
    /**
     * It sets every node to its equilibrium potential.
     * TODO: find a better guess...
     */
    void guess_equilibrium(NonlinearImplicitSystem& poisson) const;

    //! Assembles the residual vector or the jacobian matrix
    /*!
     * Assembles the residual vector or the jacobian matrix for
     * the equation system with @c Coupling T.
     *
     * This method gets called from the underlying nonlinear solver
     * library
     */
    template <int T>
    static void assemble(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);

    
};


//
// inline member functions
// 

inline
const DD::Device&
DriftDiffusion::get_device(void) const
{
  return *_device;
}

inline
DriftDiffusion::Options&
DriftDiffusion::get_options(void)
{
  return _options;
}

inline
void
DriftDiffusion::set_options(const DriftDiffusion::Options& options)
{
  _options = options;
}

inline
void
DriftDiffusion::enable_mesh_refinement(void)
{
  _options.mesh_refinement = true;
}

inline
void
DriftDiffusion::disable_mesh_refinement(void)
{
  _options.mesh_refinement = false;
}

inline
const Scaling&
DriftDiffusion::get_scaling(void) const
{
  return _scaling;
}

inline
const std::vector<std::string>&
DriftDiffusion::get_variable_names(void) const
{
  return _variables;
}

inline
const std::vector<Number>&
DriftDiffusion::get_solution(void) const
{
  return _solution;
}

inline
unsigned int
DriftDiffusion::get_n_nonlinear_iterations(void) const
{
  return _n_nonlinear_iterations;
}

inline
double
DriftDiffusion::get_final_residual(void) const
{
  return _final_residual;
}

inline
const std::map<const BoundaryDescriptor*, double>&
DriftDiffusion::get_boundary_currents() const
{
  return _boundary_currents;
}

inline
void
DriftDiffusion::solve(void)
{
  solve(false);
}

inline
void
DriftDiffusion::assign_boundary_values(double& coeff, double& value,
        const std::vector<double>& coefficients)
{
  assert(coefficients[1] != 0.0);

  coeff = coefficients[0] / coefficients[1];
  value = coefficients[2] / coefficients[1];
}


inline
EquationSystems&
DriftDiffusion::get_equation_system(void)
{
  return *_eq_system;
}

inline
Mesh& 
DriftDiffusion::get_mesh(void) const
{
  return _device->get_mesh();
}

#endif //_DRIFTDIFFUSION_H_
