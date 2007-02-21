// $Id$

#ifndef _EXCITONTRANSPORT_H_
#define _EXCITONTRANSPORT_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "Scaling.h"
#include "Device.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"

// Libmesh includes
#include "libmesh_common.h"
#include "enum_order.h"
#include "point.h"

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
class Mesh;
class Elem;
class EquationSystems;
class NonlinearImplicitSystem;
class ExcitonProperties;

template<typename T> class NumericVector;
template<typename T> class SparseMatrix;

template<typename T> class TiberPetscNonlinearSolver;
template<typename T> class NonlinearSolver;

//! The main class to perform exciton drift-diffusion calculations
/*!
 * TODO
 * Some more details
 */
class ExcitonTransport : public SimulationInterface
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
        double nonlinear_step_tolerance;
        unsigned int nonlinear_max_iterations;
        double linear_tolerance;
        double linear_abs_tolerance;
        unsigned int linear_max_iterations;
        
        //! The line search maximum step size per grid point
        double ls_maxstep;

        //! The line search type
        int ls_type;

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


      private:
        
        friend class ExcitonTransport;
    };


    //! Destructor
    ~ExcitonTransport(void);

    //! Create an ExcitonTransport object
    static ExcitonTransport* create(void);
  
    
    /*! \copydoc SimulationInterface::create_physical_model() */
    virtual PhysicalModel*
      create_physical_model(const ModelOptions& options) const
      throw (ModelErrorException);

  

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

    /*!
     * Remember the current solution for future restart
     */
    void remember_current_solution(void);

    /*!
     * Reset to the remembered solution
     */
    void set_to_remembered_solution(void);


    /*!
     * @returns the variable names in a vector.
     */
    const std::vector<std::string>& get_variable_names(void) const;

    /*!
     * @returns the current nodal solution vector.
     *
     * The vector is ordered as
     *   [ var1\@node1 ... varn\@node1 var1\@node2 ... varn\@node2 ... ]
     * where the ordering \p var1 ... \p var2 is the same as in the
     * method \p get_variable_names()
     */
    const std::vector<Number>& get_solution(void) const;

    /*!
     * @returns the number of nonlinear iterations needed for the solution
     */
    unsigned int get_n_nonlinear_iterations(void) const;

    /*!
     * @returns the final residual norm of the solution
     */
    double get_final_residual(void) const;

    //! Get the exciton electro-chemical potential in a given point
    /*!
     * If \c elem is not in the list of active elements for this simulation,
     * a parent or children which contains \p will be looked for.
     * 
     * \param elem the pointer to the element
     * \param p the point in which to calculate the potentials
     * \return the electro-chemical potential in \c p
     *
     * \pre
     * The element \c elem is assumed to contain the point \c p.
     *
     */
    double get_solution(const Elem* elem, const Point& p);
    


  protected:
      
    //! Constructor
    ExcitonTransport(void);
    
    //! Initialize the equation system
    /*!
     * This has to be called before any call to methods which access the
     * equation system object. So do it early.
     */
    virtual void do_init(void);

    //! Solve the drift-diffusion problem.
    /*!
     * If adaptive mesh refinement was enabled for this solver
     * run, it will be deactivated afterwards and has to be re-enabled
     * explicitly.
     */
    virtual void do_solve(void);

    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);


    /*! \copydoc SimulationInterface::build_nodal_results() */
    virtual void build_nodal_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    /*! \copydoc SimulationInterface::build_elemental_results() */
    virtual void build_elemental_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    
  private:

    // for nicer code
    typedef TiberPetscNonlinearSolver<Real> SolverClass;

    //! A static reference to \c this
    /*!
     * This is needed during matrix assembly, which is a static method.
     */
    static ExcitonTransport* _this;

    //! An internal pointer to the device
    Device* _device;


    /*!
     * The equation system for this device
     */
    EquationSystems* _eq_system;
    
    /*!
     * If @c true, the equation system needs to be rebuilt
     */
    bool _rebuild_eq_system;

    /*!
     * The @c Options to be used
     */
    Options _options;
    
    /*!
     * The scaling parameters
     */
    Scaling _scaling;

    /*!
     * The number of nonlinear iterations needed
     */
    unsigned int _n_nonlinear_iterations;

    /*!
     * The final residual norm
     */
    double _final_residual;

    //! disable the copy constructor
    ExcitonTransport(const ExcitonTransport& rhs);
    
    //! disable the copy assignment operator
    ExcitonTransport& operator=(const ExcitonTransport& rhs);
    
    //! Get the solution at the point \c p in a given element
    /*!
     * \param elem the pointer to the element
     * \param a vector containing the points in which to calculate the potential
     * \param solution a vector where the solutions will be stored
     *
     * \note
     * This implementation assumes, that \c elem is one of the active elements
     * of this simulation and that it contains the points in \c p.
     *
     */
    void get_solution_secure(const Elem* elem, const Point& p,
        double& solution);
 
    //! Get the electro-chemical potential at the points in \c p
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

    //! Do the actual assembly
    void do_assembly(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);

};


//
// inline member functions
// 

inline
ExcitonTransport*
ExcitonTransport::create(void)
{
  return new ExcitonTransport();
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
Mesh& 
ExcitonTransport::get_mesh(void) const
{
  return _device->get_mesh();
}

inline
void
ExcitonTransport::get_solution_secure(const Elem* elem, const Point& p,
        double& solution)
{
  std::vector<Point> pvec(1, p);
  std::vector<double> sol(1);

  get_solution_secure(elem, pvec, sol);

  solution = sol[0];
}


inline
void
ExcitonTransport::assemble(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{
  _this->do_assembly(x, residual, jacobian);
}



#endif //_EXCITONTRANSPORT_H_
