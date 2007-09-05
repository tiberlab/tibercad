// $Id$

#ifndef _DRIFTDIFFUSION_H_
#define _DRIFTDIFFUSION_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "DriftDiffusionDefs.h"
#include "Device.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"

// Libmesh includes
#include "libmesh_common.h"
#include "enum_order.h"
#include "enum_solver_type.h"
#include "enum_preconditioner_type.h"


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
class TiberLinearSolver;


template<typename T> class DenseMatrix;
template<typename T> class NumericVector;
template<typename T> class SparseMatrix;


//! The main class to perform standard drift-diffusion calculations
/*!
 * We solve the following system of equations:
 * \f{eqnarray*}
 * -\nabla(\epsilon_r\nabla\varphi -\mathrm{P})& = & \frac{e}{\epsilon_0}
 *     \left(-n + p + N_D^+ - N_A^-\right) \\
 * -\nabla(\mu_n n\nabla\phi_n) & = & R \\
 * -\nabla(\mu_p p\nabla\phi_p) & = & R
 * \f}
 * using appropriate models for the ionization of dopants, polarization,
cre* mobilities and recombinations.
 * 
 */
class DriftDiffusion : public SimulationInterface
{
  public:


    //! The variables that can be provided
    enum Variables
    {
      UNKNOWN = INVALID_ID,
      ELPOTENTIAL,
      QFERMIE,
      QFERMIH,
      CBANDEDGE,
      VBANDEDGE,
      BANDGAP,
      EDENSITY,
      HDENSITY,
      EMOBILITY,
      HMOBILITY
    };
      
 
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
      SG   /*!< Box integration with Scharfetter-Gummel */
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
     * points or elements.
     *
     * The values are given in Volts.
     * To get the quasi-Fermi energies one has to multiply the electro-chemical
     * potentials by -e (the potentials refer to a particle with charge -e):
     * \f[E_{F,n} = -e \phi_n\f]
     * \f[E_{F,p} = -e \phi_p\f]
     */
    class Solution
    {
      public:

        //! Constructor
        Solution(void);

        //! The electric potential
        double potential;
      
        //! The electron electro-chemical potential
        double fermi_e;
      
        //! The hole electro-chemical potential
        double fermi_h;

        //! Set all to some value
        Solution& operator=(double other);
    };


    //! A structure to hold the electric field
    class EField
    {
      public:

        //! Constructor
        EField(void);

        //! Get the field as RealGradient
        const RealGradient& get_field(void) const;

        //! Get the i-th component
        double operator()(unsigned int i) const;

        //! Set all values to \c other
        EField& operator=(double other);

      private:

        RealGradient _efield;

        friend class DriftDiffusion;
    };

    
    //! A structure to hold all components of the electrical current densities
    /*!
     * This structure is used for the queries of the current densities in
     * certain points or elements. It stores the x,y and x components of the
     * electrical current densities for electrons and holes.
     *
     * Current densities are given in \f$Acm^{-2}\f$
     */
    class Currents
    {
      public:

        //! Constructor
        Currents(void);

        //! Get the x-component of the electron current density
        double jn_x(void) const { return _jn_x; };
        
        //! Get the y-component of the electron current density
        double jn_y(void) const { return _jn_y; };
        
        //! Get the z-component of the electron current density
        double jn_z(void) const { return _jn_z; };
        
        //! Get the i-th component of the electron current density
        /*!
         * 0 -> x, 1 -> y, 2 -> z
         */
        double jn(unsigned int i = 0) const;
              

        //! Get the x-component of the hole current density
        double jp_x(void) const { return _jp_x; };
        
        //! Get the y-component of the hole current density
        double jp_y(void) const { return _jp_y; };
        
        //! Get the z-component of the hole current density
        double jp_z(void) const { return _jp_z; };
        
        //! Get the i-th component of the hole current density
        /*!
         * 0 -> x, 1 -> y, 2 -> z
         */
        double jp(unsigned int i = 0) const;
        
           
        //! Get the x-component of the total current density
        double j_x(void) const { return _jn_x + _jp_x; };
        
        //! Get the y-component of the total current density
        double j_y(void) const { return _jn_y + _jp_y; };
        
        //! Get the z-component of the total current density
        double j_z(void) const { return _jn_z + _jp_z; };
        
        //! Get the i-th component of the total current density
        /*!
         * 0 -> x, 1 -> y, 2 -> z
         */
        double j(unsigned int i = 0) const;

        // Get the absolute value of the total current density
        double j_abs(void) const;

        //! Set all to some value
        Currents& operator=(double other);
            
      private:
        
        double _jn_x;
        double _jn_y;
        double _jn_z;

        double _jp_x;
        double _jp_y;
        double _jp_z;

        friend class DriftDiffusion;
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

        //! The linear solver type
        /*!
         * The implementation of the linear solver depends on the
         * underlying library to be used.
         *
         */
        SolverType ksp_type;

        //! The preconditioner
        /*!
         * The implementation of the preconditioner depends on the
         * underlying library to be used.
         * 
         */
        PreconditionerType pc_type;

        //! The nonlinear solver type
        /*!
         * This means the implementation of the nonlinear system solver
         */
        std::string nonlinear_solver;

        //! The linear solver to be used (PETSc, PARDISO etc)
        std::string linear_solver;


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

        //! Use exact jacobian or not
        bool exact_newton;

      private:
        

        friend class DriftDiffusion;
    };
      


    //! Constructor
    DriftDiffusion(void);
    
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


    //! Get the band edges in a given element
    /*!
     * \param elem the pointer to the element
     * \param band_edges the band_edges
     * \note { \c band_edges[0] = Ec, \c band_edges[1] = Ev }
     *
     * This method returns the band edges of the material of element elem
     * without adding the electric potential.
     */
    void get_band_edges(const Elem* elem, std::vector<double>& band_edges);


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


     //!  Get electron thermoelectric power at the element centroid
    double get_electrons_thermoelectric_power(const Elem* elem);

    //!  Get hole thermoelectric power at the element centroid
    double get_holes_thermoelectric_power(const Elem* elem);

    //!  Get electron electron conducibility at the element centroid
    double get_electron_conducibility(const Elem* elem);

    //!  Get hole conducibility at the element centroid
    double get_hole_conducibility(const Elem* elem);
    



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

    
    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);

    
    /*! \copydoc SimulationInterface::get_solution_vector() */
    virtual NumericVector<double>& get_solution_vector(void);
    

    /*! \copydoc SimulationInterface::build_nodal_results() */
    virtual void build_nodal_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    /*! \copydoc SimulationInterface::build_elemental_results() */
    virtual void build_elemental_results(const std::set<std::string>& variables,
        std::vector<double>& results, std::vector<std::string>& legend);


    /*! \copydoc SimulationInterface::build_integrated_quantities() */
    virtual void build_integrated_quantities(
        const std::set<std::string>& names,
        std::vector<double>& values);


    /*! \copydoc SimulationInterface::build_integrated_quantities_description()
     */
    virtual void build_integrated_quantities_description(
        const std::set<std::string>& names,
        std::vector<std::string>& legend,
        std::vector<std::string>& description);


    /*! \copydoc SimulationInterface::do_maximum_norm_of_difference() */
    virtual double do_maximum_norm_of_difference(ID id);


    /*! \copydoc SimulationInterface::convert_variable_to_id() */
    virtual ID convert_variable_name_to_id(const std::string& variable_name);


    /* \copydoc SimulationInterface::get_solution_secure(const Elem*,
     * const std::vector<ID>&, std::vector<std::vector<double> >&)
     */
    virtual void get_solution_secure(const Elem* elem,
        const std::set<ID>& ids,
        std::vector<std::map<ID, double> >& values);


    /* \copydoc SimulationInterface::get_solution_secure(const Elem*,
     * const std::vector<Point>&, const std::vector<ID>&,
     * std::vector<std::vector<double> >&)
     */
    virtual void get_solution_secure(const Elem* elem,
        const std::vector<Point>& p,
        const std::set<ID>& ids,
        std::vector<std::map<ID, double> >& solution);


  private:
   


    // for nicer code
    typedef std::map<const Boundary*, double> ContactData;
    typedef std::map<const Node*, Boundary*> BoundaryNodeList;

    //! A static reference to \c this
    /*!
     * This is needed during matrix assembly, which is a static method.
     */
    static DriftDiffusion* _this;

    //! An internal pointer to the device
    Device* _device;


    //! A linear solver
    /*!
     * The linear solver is used to get a good guess for the electro-chemical
     * potentials
     */
    TiberLinearSolver* _linear_solver;

    /*!
     * A list of nodes with dirichlet boundary conditions
     */
    BoundaryNodeList _dirichlet_nodes;


    /*!
     * \brief A list of nodes which lie on an inner boundary between
     * a dielectric and a semiconductor
     */
    std::set<const Node*> _dielectric_boundary_nodes;


    /*!
     * If @c true, the equation system needs to be rebuilt
     */
    bool _rebuild_eq_system;


    /*!
     * The @c Options to be used
     */
    Options _options;
    

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
  
    
    
    //! Rebuild the equation system if needed
    void rebuild_equation_system(void);

    
    /*!
     * \brief Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(Scaling::ScalingType type = Scaling::UNITS);

    
    
    //! Fills the dirichlet nodes data structure.
    void find_dirichlet_nodes(void);


    //! Find nodes on boundary between dielectric/semiconductor
    void find_dielectric_boundary_nodes(void);
    

    //! Tells if node lies on an inner dielectric/semiconductor boundary
    bool is_dielectric_boundary_node(const Node* node) const;

    
    //! Reset solver environment.
    /*!
     * Deletes only the \p EquationSystems object, without
     * touching simulation voltages and solutions.
     */
    void reset_solver(void);


    //! Sets the Dirichlet type boundary conditions
    void set_dirichlet_bc(void);
    

    //! Cleanup solver environment.
    /*!
     * Deletes the \p EquationSystems object, the simulation voltages
     * vector and the solution and variable vectors.
     */
    void cleanup_solver(void);


    //! Solve using Newton method
    void solve_newton(void);

    //! Solve using an iterative Gummel scheme
    void solve_gummel(void);


    //! Do a Newton type iteration
    void do_newton(void);

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

    //! Get a solution at the points in \c p in a given element
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
    template <typename T>
    void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
        std::vector<T>& solution);



    //! Get the band edges of a given element
    /*!
     * \note
     * This implementation assumes, that \c elem is one of the active elements
     * of this simulation and that it contains the points in \c p.
     */
    void get_bands_secure(const Elem* elem, std::vector<double>& band_edges);



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
    template <int T>
    void do_assembly(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);
    void assemble_linear_electrons(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);
    void solve_linear(void);


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



inline
bool
DriftDiffusion::is_dielectric_boundary_node(const Node* node) const
{
  bool result = false;
  if (_dielectric_boundary_nodes.find(node) != _dielectric_boundary_nodes.end())
    result = true;

  return result;
}



inline
DriftDiffusion::Solution::Solution(void)
  : potential(0.0),
    fermi_e(0.0),
    fermi_h(0.0)
{
}



inline
DriftDiffusion::Solution&
DriftDiffusion::Solution::operator=(double other)
{
  potential = other;
  fermi_e = other;
  fermi_h = other;
  return *this;
}



inline
DriftDiffusion::Currents::Currents(void)
  : _jn_x(0.0),
    _jn_y(0.0),
    _jn_z(0.0),
    _jp_x(0.0),
    _jp_y(0.0),
    _jp_z(0.0)
{
}



inline
double
DriftDiffusion::Currents::jn(unsigned int i) const
{
  switch (i)
  {
    case 1:
      return _jn_y;
      break;
    case 2:
      return _jn_z;
      break;
    default:
      return _jn_x;
  }
}



inline
double
DriftDiffusion::Currents::jp(unsigned int i) const
{
  switch (i)
  {
    case 1:
      return _jp_y;
      break;
    case 2:
      return _jp_z;
      break;
    default:
      return _jp_x;
  }
}



inline
double
DriftDiffusion::Currents::j(unsigned int i) const
{
  switch (i)
  {
    case 1:
      return j_y();
      break;
    case 2:
      return j_z();
      break;
    default:
      return j_x();
  }
}


inline
double
DriftDiffusion::Currents::j_abs(void) const
{
  return std::sqrt(j_x() * j_x() + j_y() * j_y() + j_z() * j_z());
}



inline
DriftDiffusion::Currents&
DriftDiffusion::Currents::operator=(double other)
{
  
  _jn_x = other;
  _jn_y = other;
  _jn_z = other;
  _jp_x = other;
  _jp_y = other;
  _jp_z = other;
  return *this;
}





inline
DriftDiffusion::EField::EField(void)
  : _efield(0.0)
{
}


inline
const RealGradient&
DriftDiffusion::EField::get_field(void) const
{
  return _efield;
}


inline
DriftDiffusion::EField&
DriftDiffusion::EField::operator=(double other)
{
  _efield = other;

  return *this;
}


inline
double
DriftDiffusion::EField::operator()(unsigned int i) const
{
  return _efield(i);
}




#endif //_DRIFTDIFFUSION_H_
