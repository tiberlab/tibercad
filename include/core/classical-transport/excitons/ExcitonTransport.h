// $Id$

#ifndef _EXCITONTRANSPORT_H_
#define _EXCITONTRANSPORT_H_

#include "SimulationOptions.h"
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
#include <map>

// forward declarations
class DD::Device;
class Mesh;
class Elem;
class Node;
class EquationSystems;
class NonlinearImplicitSystem;
class DriftDiffusion;
class ExcitonProperties;

template<typename T> class NumericVector;
template<typename T> class SparseMatrix;

template<typename T> class TiberPetscNonlinearSolver;
template<typename T> class NonlinearSolver;

//! The main class to perform standard drift-diffusion calculations
/*!
 * TODO
 * Some more details
 */
class ExcitonTransport
{
  public:
 
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

        friend class ExcitonTransport;
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
         * The order of gauss integration
         */
        libMeshEnums::Order integration_order;

        /**
         * The nonlinear/linear (PETSc-)solver parameters
         */
        SolverParameters solver_params;

        /**
         * The units in which the mesh object is given
         */
        double mesh_units;

        //! Use a local (nodal based) scaling for continuity equations
        /*!
         * This probably should be the default and should replace the global
         * scaling factors for the continuity equations
         */
        bool local_scaling;

      private:
        
        friend class ExcitonTransport;
    };

      
    ExcitonTransport(DD::Device* device);

    ExcitonTransport(DD::Device* device, ExcitonTransport::Options& params);

    ~ExcitonTransport(void);
    
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

    //! Set the exciton model
    void set_exciton_model(ExcitonProperties* exciton_model)
      { _exciton_model = exciton_model; };

    //! Get a pointer to the exciton model
    ExcitonProperties* get_exciton_model(void)
      { return _exciton_model; };

    /**
     * @returns a reference to the simulation options
     */
    Options& get_options(void);

    //! Initialise the equation system
    void init(void);
    
    //! Set the DriftDiffusion object
    void set_driftdiffusion(DriftDiffusion* dd)
      { _dd_object = dd; };

    //! Get the DriftDiffusion object
    DriftDiffusion* get_driftdiffusion(void)
      { return _dd_object; };

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
    
    //! Set an initial guess
    void set_initial_guess(double guess);

    
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
    //void build_scaling(void);

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
     * fill @c densities with electron and hole densities
     *
     * TODO add other values
     */
    void build_densities(std::vector<double>& densities,
        std::vector<std::string>& names);

    //! Fill a vector with the electric field data
    void build_current_density(std::vector<double>& current,
        std::vector<std::string>& names);

  private:

    // for nicer code
    typedef TiberPetscNonlinearSolver<Real> SolverClass;

    //! A static reference to \c this
    /*!
     * This is needed during matrix assembly, which is a static method.
     */
    static ExcitonTransport* _this;

    /**
     * The device to be solved by this ExcitonTransport object
     */
    DD::Device* _device;

    DriftDiffusion* _dd_object;

    ExcitonProperties* _exciton_model;

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
     * The number of nonlinear iterations needed
     */
    unsigned int _n_nonlinear_iterations;

    /**
     * The final residual norm
     */
    double _final_residual;

    //! disable the copy constructor and assignment operator
    ExcitonTransport(const ExcitonTransport& rhs);
    ExcitonTransport& operator=(const ExcitonTransport& rhs);

    /**
     * Set the options for the PETSc solver as given in @c SolverParameters
     */
    void set_solver_params(NonlinearSolver<Number>& solver);

    /**
     * Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(void);

    //! Reset solver environment.
    /*! 
     * Deletes the \p EquationSystems object, the simulation voltages
     * vector and the solution and variable vectors.
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

    //! Assembles the residual vector or the jacobian matrix
    /*!
     * Assembles the residual vector or the jacobian matrix for
     * the equation system with @c Coupling T.
     *
     * This method gets called from the underlying nonlinear solver
     * library
     */
    static void assemble(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);

};


//
// inline member functions
// 

inline
const DD::Device&
ExcitonTransport::get_device(void) const
{
  return *_device;
}

inline
ExcitonTransport::Options&
ExcitonTransport::get_options(void)
{
  return _options;
}

inline
void
ExcitonTransport::set_options(const ExcitonTransport::Options& options)
{
  _options = options;
}

inline
void
ExcitonTransport::enable_mesh_refinement(void)
{
  _options.mesh_refinement = true;
}

inline
void
ExcitonTransport::disable_mesh_refinement(void)
{
  _options.mesh_refinement = false;
}

inline
const Scaling&
ExcitonTransport::get_scaling(void) const
{
  return _scaling;
}

inline
unsigned int
ExcitonTransport::get_n_nonlinear_iterations(void) const
{
  return _n_nonlinear_iterations;
}

inline
double
ExcitonTransport::get_final_residual(void) const
{
  return _final_residual;
}

inline
EquationSystems&
ExcitonTransport::get_equation_system(void)
{
  return *_eq_system;
}

inline
Mesh& 
ExcitonTransport::get_mesh(void) const
{
  return _device->get_mesh();
}

#endif //_EXCITONTRANSPORT_H_
