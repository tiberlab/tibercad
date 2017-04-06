
#ifndef _MASTEREQUATION_H_
#define _MASTEREQUATION_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "MasterEquationDefs.h"
#include "TiberNonlinearSystem.h"

#include "TiberModelObject.h"
#include "Boundary.h"


// Libmesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/linear_implicit_system.h"


// C++ includes
#include <vector>
#include <set>
#include <map>

class Database;
class PhysicalObject;


// forward declarations
namespace libMesh
{

class MeshBase;
//class Device;
//class Boundary;
//class MeshBase;
//class Elem;
//class Point;
class Node;
class EquationSystems;



template<typename T> class DenseMatrix;
template<typename T> class NumericVector;
template<typename T> class SparseMatrix;

}

class TiberLinearSolver;


class TBDLLOCAL MasterEquation : public SimulationInterface
{

  public:

    enum Solutions
    {
      eFermiLevel,        /*!< electron density */
      hFermiLevel,        /*!< hole density */
      LUMOLevel,         /*!< conduction band lowest edge */
      HOMOLevel,          /*!< valence band highest edge */
      eOccProbability,         /*!< conduction band lowest edge */
      hOccProbability,
      //Potential,        /*!< the potential */
      //ElField,            /*!< the field (negative gradient of potential) */
      //Displacement,     /*!< the electric displacement */
      //Polarization,
      ChargeDensity     /*!< the source (charge density) */
    };

    //! The solver methods that can be used
    enum SolverMethod
    {
      NEWTON
    };

    /**
     * This class defines various parameters that control a
     * Master Equations calculation
     */

    class Options
    {
      public:

        Options(void);

        Options(const Options& rhs);

        Options& operator=(const Options& rhs);

        int coupling;

        SolverMethod solver_method;

        bool _xmonitor = true;


      private:

        friend class MasterEquation;

    };


    // Destructor
    virtual ~MasterEquation(void);


    //! Create an DriftDiffusion object
    static MasterEquation* create(const ModelOptions& options);


/*
    //! We need to create a physical model
    virtual PhysicalModel*
      create_bulk_model(const ModelOptions& options,
        const Material* mat) const;
*/

    /**
     * @returns a reference to the simulation options
     */
    Options& get_options(void);


    //! Get the mesh
    /*!
     * \return a constant reference to the simulation mesh
     */
    //MeshBase& get_mesh(void) const;


    //! Get the number of nonlinear iterations needed for the solution
    unsigned int get_n_nonlinear_iterations(void) const;


    //! Get the final residual norm of the solution
    double get_final_residual(void) const;


	const Database& get_database(void);



  protected:

    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    MasterEquation(const ModelOptions& options);



    //! The initialization
    virtual void do_init(void);


    //! Read database
	virtual void read_database(void);



    //! Parse the options from the input file
    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the Master Equations
    virtual void do_solve(void);


    //! Print some useful information
    virtual void do_print_info(void);



    //! We have to provide somehow our solution variables // DA MODIF
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& solutions,
        const std::vector<Point>& points);



    //! We override this to not produce too many files ...
    virtual void plot_globaldata(void) {};

  private:

    void cleanup_solver(void);

    void reset_solver(void);

    //! A static pointer to this
    static MasterEquation* _this;

    Device* _device;
	

	//! last added
	SimulationEnvironment* _env;
	
	
	//! last added
	std::map<const Boundary*, unsigned int> _bd_num;
	
	
	
	const PhysicalObject* _owner;


    //! Do a Newton type iteration
    void do_newton(void);
	
	
	


	

    /*!
     * If @c true, the equation system needs to be rebuilt
     */
    bool _rebuild_eq_system;


    unsigned int _n_nonlinear_iterations;

    /**
     * The final residual norm
     */
    double _final_residual;



    /*!
     * The @c Options to be used
     */
    Options _options;


    /**
     * @returns a reference to the simulation options
     */
    Options& get_my_options(void);


    //! The assembly function
    //static void assemble(EquationSystems& es, const std::string& system_name);

   static void assemble(const libMesh::NumericVector<Number>& x,
       libMesh::NumericVector<Number>* residual,
         libMesh::SparseMatrix<Number>* jacobian,
           libMesh::NonlinearImplicitSystem& sys);


   static void transformation(libMesh::NumericVector<Number>& u,
       libMesh::NumericVector<Number>& T,
         libMesh::NumericVector<Number>& TX, bool transf,
           libMesh::NonlinearImplicitSystem&);


    //! The real assembly function
    //void do_assembly(EquationSystems& es, const std::string& system_name);



    template <int T>
    void do_assembly(const libMesh::NumericVector<Number>& x, libMesh::NumericVector<Number>* residual,
        libMesh::SparseMatrix<Number>* jacobian);


    void do_transformation(libMesh::NumericVector<Number>& u,
        libMesh::NumericVector<Number>& T,
          libMesh::NumericVector<Number>& TX, bool transf);


    //! Rebuild the equation system if needed
    void rebuild_equation_system(void);



    //TiberLinearSystem *_sys_EcEv;

};


inline
MasterEquation*
MasterEquation::create(const ModelOptions& options)
{
  return new MasterEquation(options);
}


inline
MasterEquation::Options&
MasterEquation::get_my_options(void)
{
  return _options;
}


inline
MasterEquation::Options&
MasterEquation::get_options(void)
{
  return _options;
}


inline
unsigned int
MasterEquation::get_n_nonlinear_iterations(void) const
{
  return _n_nonlinear_iterations;
}

inline
double
MasterEquation::get_final_residual(void) const
{
  return _final_residual;
}

/**
inline
 const PhysicalObject*
 MasterEquation::get_owner(void) const
 {
   return _owner;
 }
*/



#endif // _CONTINUOUSKINETIC_H_
