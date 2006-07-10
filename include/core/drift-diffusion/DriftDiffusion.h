// $Id$

#ifndef _DRIFTDIFFUSION_H_
#define _DRIFTDIFFUSION_H_

#include "SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "Scaling.h"
#include "DDevice.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"

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
#include <set>
#include <map>

// forward declarations
class DD::Device;
class ElectricalContact;
class Mesh;
class Elem;
class Point;
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

    //! A structure to hold the three potentials
    /*!
     * This structure is used for the queries of the solution in certain
     * points or elements
     */
    struct Solution
    {
      //! The electric potential
      double potential;
      
      //! The electron electro-chemical potential
      double fermi_e;
      
      //! The hole electro-chemical potential
      double fermi_h;
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

        //! Include artificial drift in continuity equations
        /*!
         * This can be helpful in materials of high bandgap, but it is
         * a somewhat dirty trick.
         */
        bool artificial_drift;

        //! Use a local (nodal based) scaling for continuity equations
        /*!
         * This probably should be the default and should replace the global
         * scaling factors for the continuity equations
         */
        bool local_scaling;

        //! linearize continuity equations
        bool linearize_continuity_eq;


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

    void set_equation_systems(EquationSystems* eq_systems)
      { _eq_system = eq_systems; }

    //! Initialize the equation system
    /*!
     * This has to be called before any call to methods which access the
     * equation system object. So do it early.
     */
    void init(void);

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

    //! Calculate the nodal scaling factors for the continuity equations
    /*!
     * This method currently calculates the densities on each node as scaling
     * factors for the electron and hole continuity equations.
     * A better way would be to use some potential based values which can be
     * calculated from nodal potential values during matrix assembly.
     */
    void build_scaling(void);

    /**
     * Set the simulation voltage for the boundary named \p boundary.
     */
    // TODO throw exception if boundary non-existing
    void set_simulation_voltage(const std::string& boundary,
        double voltage);

    //! Remember the current solution for future restart
    void remember_current_solution(void);

    //! Reset to the remembered solution
    void set_to_remembered_solution(void);

    //! Set the electron quasi Fermi level to \c Ef_n
    void set_electron_fermi_level(double Ef_n);

    //! Set the hole quasi Fermi level to \c Ef_p
    void set_hole_fermi_level(double Ef_p);

    //! Solve the drift-diffusion problem.
    /*!
     * If adaptive mesh refinement was enabled for this solver
     * run, it will be deactivated afterwards and has to be re-enabled
     * explicitly.
     */
    void solve(void);

    /**
     * @returns the current nodal solution vector.
     *
     * The vector is ordered as
     *   [ var1\@node1 ... varn\@node1 var1\@node2 ... varn\@node2 ... ]
     * where the ordering \p var1 ... \p var2 is the same as in the
     * method \p get_variable_names()
     */
    const std::vector<double>& get_solution(void) const;

    //! Get the solution on the nodes of a given element
    /*!
     * \param elem the pointer to the element
     * \param[inout] solution the vector where the solution will be stored
     */
    void get_solution(const Elem* elem, std::vector<Solution>& solution);

    //! Get the solution at the point p in a given element
    /*!
     * If \c elem is not in the list of active elements for this simulation,
     * a parent or children which contains \p will be looked for.
     * 
     * \param elem the pointer to the element
     * \param p the point in which to calculate the potentials
     * \param[inout] solution the structure where the solution will be stored
     *
     * \pre
     * The element \c elem is assumed to contain the point \c p.
     *
     */
    template <typename T>
    void get_solution(const Elem* elem, const Point& p, T& solution);

    //! Get the solution at the points p in a given element
    /*!
     * If \c elem is not in the list of active elements for this simulation,
     * a parent or children which contains the points in \c p
     * will be looked for.
     * 
     * \param elem the pointer to the element
     * \param p vector with the points in which to calculate the potentials
     * \param[inout] solution the vector where the solutions
     * will be stored
     *
     */
    template <typename T>
    void get_solution(const Elem* elem, const std::vector<Point>& p,
        std::vector<T>& solution);


    //! Get the electric potential at a given point in a given element
    /*!
     * \param elem the pointer to the element
     * \param p  vector with the points in which to calculate the potential
     * \param[inout] potential the vector where the potential will be stored
     */
    void get_electric_potential(const Elem* elem, const std::vector<Point>& p,
        std::vector<double>& potential);


    //! Get the electric potential at a given point in a given element
    /*!
     * \param elem the pointer to the element
     * \param p the point in which the potential should be calculated
     * \return the potential in point \c p
     */
    double get_electric_potential(const Elem* elem, const Point& p);

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
    const std::map<const ElectricalContact*, double>&
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

    //! Fill a vector with the electric field data
    void build_current_density(std::vector<double>& current,
        std::vector<std::string>& names);

    //! Fill a vector with the electric field data
    void build_electric_field(std::vector<double>& field,
        std::vector<std::string>& names);

    //! Fill a vector with the band edge data
    void build_band_edges(std::vector<double>& band_edges,
        std::vector<std::string>& names);

    //! Fill a vector with the elemental band edge data
    void build_elem_band_edges(std::vector<double>& band_edges,
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
    typedef std::map<const ElectricalContact*, double> ContactData;
    typedef std::map<const Node*, ElectricalContact*> BoundaryNodeList;
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

    //! The elements that were used in the last simulation
    std::set<const Elem*> _element_list;

    //! Update the element list
    void update_element_list(void);

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
    std::vector<double> _solution;

    /**
     * The number of nonlinear iterations needed
     */
    unsigned int _n_nonlinear_iterations;

    /**
     * The final residual norm
     */
    double _final_residual;

    //! disable the copy constructor and assignment operator
    DriftDiffusion(const DriftDiffusion& rhs);
    DriftDiffusion& operator=(const DriftDiffusion& rhs);

    //! Do a number of Gummel iterations
    /*!
     * \return the final residual
     * \param the maximum number of iterations
     */
    double do_gummel_iterations(int max_it)
      throw (PetscRuntimeError, KSPDivergedError, SNESDivergedError);

    /**
     * Set the options for the PETSc solver as given in @c SolverParameters
     */
    void set_solver_params(NonlinearSolver<Number>& solver,
        CalculationType calc_type = NONEQUILIBRIUM);

    /**
     * Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(Scaling::ScalingType type = Scaling::UNITS);

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


    //! Get the equation system
    EquationSystems& get_equation_system(void);

    //! Solve using Newton method
    void solve_newton(void) throw (PetscRuntimeError);

    //! Solve using an iterative Gummel scheme
    void solve_gummel(void) throw (PetscRuntimeError);

    //! Calculate terminal currents
    void calculate_currents(void);

    //! Get the solution at the point \c p in a given element
    /*!
     * Depending on the template argument, get the electric potential
     * or all three potentials.
     *
     * \param elem the pointer to the element
     * \param p the point in which to calculate the potentials
     * \param solution a reference to the structure where the solution will be
     * put into
     *
     * \note
     * This implementation assumes, that \c elem is one of the active elements
     * of this simulation and that it contains \c p.
     *
     */
    template <typename T>
    void get_solution_secure(const Elem* elem, const Point& p,
        T& solution);

    //! Get the solution at the points in \c p in a given element
    /*!
     * \param elem the pointer to the element
     * \param a vector containing the points in which to calculate the potentials
     * \param solution a vector where the solutions will be stored
     *
     * \note
     * This implementation assumes, that \c elem is one of the active elements
     * of this simulation and that it contains the points in \c p.
     *
     */
    void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
        std::vector<Solution>& solution);

    //! Get the electric potential at the points in \c p in a given element
    /*!
     * \param elem the pointer to the element
     * \param a vector containing the points in which to calculate the potentials
     * \param solution a vector where the potential will be stored
     *
     * \note
     * This implementation assumes, that \c elem is one of the active elements
     * of this simulation and that it contains the points in \c p.
     *
     */
    void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
        std::vector<double>& solution);



    void build_solution_vector(std::vector<double>& vector);

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
     */
    void guess_equilibrium(void);

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
const std::vector<double>&
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
const std::map<const ElectricalContact*, double>&
DriftDiffusion::get_boundary_currents() const
{
  return _boundary_currents;
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

template <typename T>
inline
void
DriftDiffusion::get_solution_secure(const Elem* elem, const Point& p,
        T& solution)
{
  std::vector<Point> pvec(1, p);
  std::vector<T> sol(1);

  get_solution_secure(elem, pvec, sol);

  solution = sol[0];
}

inline
double
DriftDiffusion::get_electric_potential(const Elem* elem, const Point& p)
{
  double pot;
  get_solution(elem, p, pot);

  return pot;
}

inline
void
DriftDiffusion::get_electric_potential(const Elem* elem,
    const std::vector<Point>& p, std::vector<double>& potential)
{
  get_solution(elem, p, potential);
}




#endif //_DRIFTDIFFUSION_H_
