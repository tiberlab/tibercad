// $Id$

#ifndef _DRIFTDIFFUSION_H_
#define _DRIFTDIFFUSION_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "Device.h"
#include "Scaling.h"
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
//class Device;
class Boundary;
class Mesh;
class Elem;
class Point;
class Node;
class EquationSystems;
class NonlinearImplicitSystem;

template<typename T> class DenseMatrix;
template<typename T> class NumericVector;
template<typename T> class SparseMatrix;

template<typename T> class TiberPetscNonlinearSolver;
template<typename T> class NonlinearSolver;

//! The main class to perform standard drift-diffusion calculations
/*!
 * TODO
 * Some more details
 */
class DriftDiffusion : public SimulationInterface
{
  public:
 
    //! The solver methods that can be used
    enum SolverMethod
    {
      NEWTON,
      GUMMEL
    };

    //! The implemented discretization schemes
    enum DiscretizationScheme
    {
      FEM, /*!< Standard Finite Elements */
      BOX, /*!< Box integration method */
      FEMVARIANT /*!< Finite Elements on slightly different equations */
    };

    //! How to calculate currents
    enum CurrentCalculation
    {
      /*!
       * Use Ramo-Shockley test functions
       * (cf. calculate_current_rstf() )
       */
      RSTF,
      /*!
       * Do a naive surface integration
       * (cf. calculate_current_surfint() )
       */
      SURFINT
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

        //! The type of scaling to be applied
        /*!
         * Can be one of \c NONE, \c UNITS or \c DEMARI
         */
        Scaling::ScalingType scaling_type;

        //! The type of coupling to be used
        /*!
         * We safe it as int, because we want to assign values by
         * using logic operators, e.g. \code ELECTRONS | POISSON \endcode
         */
        int coupling;

        //! The discretization method
        DiscretizationScheme scheme;

        //! How to calculate currents
        CurrentCalculation current_calculation; 


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
      

    
    //! Destructor
    virtual ~DriftDiffusion(void);

    
    //! Create an DriftDiffusion object
    static DriftDiffusion* create(void);

    
    /*! \copydoc SimulationInterface::create_physical_model() */
    virtual PhysicalModel*
      create_physical_model(const ModelOptions& options) const
      throw (ModelErrorException);


    /*! \copydoc SimulationInterface::create_boundary_model() */
    virtual BoundaryProperties*
      create_boundary_model(const ModelOptions& options) const
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
    
    
    /**
     * @returns a constant reference to the current scaling parameters
     */
    const Scaling& get_scaling(void) const;


    //! Set the electron quasi Fermi level to \c Ef_n
    /*!
     * \c Ef_n has to be given as electron energy in units of eV
     *
     * \note {Put some figure her for illustration. }
     */
    void set_electron_fermi_level(double Ef_n);

    
    //! Set the hole quasi Fermi level to \c Ef_p
    /*!
     * \c Ef_p has to be given as electron energy in units of eV
     *
     * \note {Put some figure her for illustration. }
     */
    void set_hole_fermi_level(double Ef_p);

    
    //! Set the electric potential everywhere to \c pot
    void set_electric_potential(double pot);

    
    //! Makes a first guess of the equilibrium potential
    /**
     * It sets every node to its equilibrium potential.
     */
    void guess_equilibrium(void);

    
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

    
    //! Get the number of nonlinear iterations needed for the solution
    unsigned int get_n_nonlinear_iterations(void) const;

    
    //! Get the final residual norm of the solution
    double get_final_residual(void) const;

    
    //! Get the boundary currents indexed by boundary descriptor
    const std::map<const Boundary*, double>&
      get_boundary_currents(void) const;


    //! Build a vector with all densities and recombination rates
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



  protected:

    //! Constructor
    DriftDiffusion(void);


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



  private:

    // for nicer code
    typedef std::map<const Boundary*, double> ContactData;
    typedef std::map<const Node*, Boundary*> BoundaryNodeList;
    typedef TiberPetscNonlinearSolver<Real> SolverClass;

    //! A static reference to \c this
    /*!
     * This is needed during matrix assembly, which is a static method.
     */
    static DriftDiffusion* _this;

    //! An internal pointer to the device
    Device* _device;

    /*!
     * A list of nodes with dirichlet boundary conditions
     */
    BoundaryNodeList _dirichlet_nodes;


    /*!
     * If @c true, the equation system needs to be rebuilt
     */
    bool _rebuild_eq_system;

    //! Decides if only equilibrium has to be solved
    bool _solve_equilibrium;


    /*!
     * The @c Options to be used
     */
    Options _options;
    
    //! The scaling parameters
    Scaling _scaling;

    //! The boundary currents
    /*!
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


    //! disable the copy constructor
    DriftDiffusion(const DriftDiffusion& rhs);
    
    //! disable the copy assignment operator
    DriftDiffusion& operator=(const DriftDiffusion& rhs);
    

    //! Do a number of Gummel iterations
    /*!
     * \return the final residual
     * \param the maximum number of iterations
     */
    double do_gummel_iterations(int max_it)
      throw (PetscRuntimeError, KSPDivergedError, SNESDivergedError);

    //! Parse the options which will not change between calls to solve()
    void parse_const_options(void);
    
    /**
     * Set the options for the PETSc solver as given in @c SolverParameters
     */
    void set_solver_params(NonlinearSolver<Number>& solver);

    //! Rebuild the equation system if needed
    void rebuild_equation_system(void);

    /**
     * Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(Scaling::ScalingType type = Scaling::UNITS);

    /**
     * Fills the dirichlet nodes data structure.
     */
    void find_dirichlet_nodes(void);

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


    //! Solve using Newton method
    void solve_newton(void);

    //! Solve using an iterative Gummel scheme
    void solve_gummel(void) throw (PetscRuntimeError);


    //! Solve the equilibrium
    void solve_equilibrium(void);


    //! Calculate terminal currents
    /*!
     * Integrates numerically over the boundary elements.
     * The current on the l'th contact is then:
     *
     * \f[I_l = -e \int_{\Gamma_l}\left(\mu_n n \nabla\phi_n +
     * \mu_p p \nabla\phi_p \right) \mathrm{d}\Gamma \f]
     */
    void calculate_currents_surfint(void);

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



    //! Assemble the residual vector or the jacobian matrix
    /*!
     * This method gets called from the underlying nonlinear solver
     * library. It calls the real assembly routines.
     */
    static void assemble(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);


    //! Assembles the residual vector or the jacobian matrix
    /*!
     * Assembles the residual vector or the jacobian matrix for
     * the equation system with @c Coupling T.
     *
     * This implementation uses standard FEM.
     */
    template <int T>
    void do_assembly(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);
 

    //! Assembles the residual vector or the jacobian matrix
    /*!
     * Assembles the residual vector or the jacobian matrix for
     * the equation system with @c Coupling T.
     *
     * This implementation uses standard FEM, but with the continuity
     * equations written in a different form.
     */
    template <int T>
    void do_assembly_new(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);
    
   
    //! Assembles the residual vector or the jacobian matrix for 1D
    /*!
     * Assembles the residual vector or the jacobian matrix for
     * the equation system with @c Coupling T.
     *
     * This implementation is for 1D only and implements the 
     * Box Integration Method.
     */
    template <int T>
    void do_assembly1D(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);

    void make_Mmatrix(DenseMatrix<Number>& m);

};


//
// inline member functions
// 

inline
DriftDiffusion*
DriftDiffusion::create(void)
{
  return new DriftDiffusion();
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
const std::map<const Boundary*, double>&
DriftDiffusion::get_boundary_currents() const
{
  return _boundary_currents;
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
