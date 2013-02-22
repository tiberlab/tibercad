// $Id$

#ifndef _EXCITONTRANSPORT_H_
#define _EXCITONTRANSPORT_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "Device.h"
#include "PetscRuntimeError.h"
#include "KSPDivergedError.h"
#include "SNESDivergedError.h"
#include "ExcPolProps.h"
#include "TiberNonlinearScalarSolver.h"
// Libmesh includes
#include "libmesh_common.h"
#include "enum_order.h"
#include "point.h"


// C++ includes
#include <vector>
#include <map>

// forward declarations
class MeshBase;
class Elem;
class EquationSystems;
class ExcitonProperties;

template<typename T> class NumericVector;
template<typename T> class SparseMatrix;

template<typename T> class NonlinearSolver;

//! The main class to perform exciton drift-diffusion calculations
/*!
 * We solve the following equation:
 * \f[-\nabla(\mu_x x\nabla\phi_x)  =  R\f]
 *
 * The get_solution() methods can provide the following variables:
 * \li \c xDensity exciton density (cm^-3)
 * \li \c xCond exciton conductivity (S/cm)
 * \li \c xMob exciton mobility (cm^2/Vs)
 * \li \c chemPot exciton chemical potential (eV)
 * \li \c J modulus of total exciton flux (cm^-2)
 * \li \c J_x exciton flux density, x-component
 * \li \c J_y exciton flux density, y-component
 * \li \c J_z exciton flux density, z-component
 *
 */
class ExcitonTransport : public SimulationInterface
{
  public:

    //! The variables that can be provided
    enum Variables
    {
      UNKNOWN = 0,
      CHEMPOT,          /*!< chemical potential */
      XDENSITY,         /*!< density */
      XMOBILITY,        /*!< mobility */
      XSIGMA,           /*!< conductivity */
      J,                /*!< total flux, modulus */
      RDISS,             /*!< dissociation rate */
      RRAD,             /*!< dissociation rate */
      RNONRAD,             /*!< dissociation rate */
      EXCEXCPOLARITON,             /*!< dissociation rate */
      EXCPHONONPOLARITON,
      NETRECOMB,             /*!< dissociation rate */
      NPOLARITON,        /* Polaritons number */
      XGEN,              /* Exciton generation */
      RADPOWER,        /*!<Radiative power density  */
      POLRADPOWER,     /* polariton radiation power */
      MOTTDENSITY           /* max.density/maximum allowed density*/
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


      private:
        
        friend class ExcitonTransport;
    };

    

    //! Constructor
    ExcitonTransport(const ModelOptions& options);

    //! Destructor
    ~ExcitonTransport(void);

    //! Create an ExcitonTransport object
    static ExcitonTransport* create(const ModelOptions& options);
  
    
    /*! \copydoc SimulationInterface::create_bulk_model() */
    virtual PhysicalModel*
      create_bulk_model(const ModelOptions& options,
          const Material* mat) const;

  

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
    MeshBase& get_mesh(void) const;
    
    //! Set an initial guess
    void set_initial_guess(double guess);

    
    //! Calculate the nodal scaling factors for the continuity equations
    /*!
     * This method currently calculates the densities on each node as scaling
     * factors for the electron and hole continuity equations.
     * A better way would be to use some potential based values which can be
     * calculated from nodal potential values during matrix assembly.
     */
    //void build_scaling(void);


    /*!
     * @returns the number of nonlinear iterations needed for the solution
     */
    unsigned int get_n_nonlinear_iterations(void) const;

    /*!
     * @returns the final residual norm of the solution
     */
    double get_final_residual(void) const;


    ExcPolProps* get_excpol_props();

  protected:
    
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

    /*! \copydoc SimulationInterface::do_maximum_norm_of_difference() */
    virtual double do_maximum_norm_of_difference(ID id);


    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

    virtual void do_setup_solution_variables(void);

    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);

    //! We override this to not produce too many files ...
    virtual void plot_globaldata(void) {};
    
    virtual void calc_excpol_integral(const NumericVector<Number>& x);
  private:
    TiberNonlinearScalarSolver scalarSolver;

    bool polaritons;

    ExcPolProps excpolprops;

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
    
   
    /**
     * Computes the scaling parameters according to the
     * scaling type \p type
     */
    void compute_scaling(void);


    //! Calculate the local scaling values
    void build_local_scaling(void);


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

    static void assemble_DRSystem(const double& x,
        double* rhs,
        double* jacobian);

    //! Do the actual assembly
    void do_assembly(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);

    void do_assembly_DRSystem(const double& x,
        double* rhs,
        double* jacobian);
};


//
// inline member functions
// 

inline
ExcitonTransport*
ExcitonTransport::create(const ModelOptions& options)
{
  return new ExcitonTransport(options);
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
MeshBase& 
ExcitonTransport::get_mesh(void) const
{
  return _device->get_mesh();
}



inline
void
ExcitonTransport::assemble(const NumericVector<Number>& x,
    NumericVector<Number>* residual,
    SparseMatrix<Number>* jacobian)
{
  _this->do_assembly(x, residual, jacobian);
}

inline
void
ExcitonTransport::assemble_DRSystem(const double& x,
    double* rhs,
    double* jacobian)
{
  _this->do_assembly_DRSystem(x, rhs, jacobian);
}

inline ExcPolProps* ExcitonTransport::get_excpol_props() {
  return &excpolprops;
}



#endif //_EXCITONTRANSPORT_H_
