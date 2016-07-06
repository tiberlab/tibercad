

#ifndef _MASTEREQUATIONS_H_
#define _MASTEREQUATIONS_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "MasterEquationsDefs.h"
#include "MasterEquationsModelInterface.h"


// Libmesh includes
#include "libmesh_common.h"
#include "enum_order.h"
#include "enum_quadrature_type.h"
#include "linear_implicit_system.h"


// C++ includes
#include <vector>
#include <set>
#include <map>


// forward declarations
class Device;
class Boundary;
class MeshBase;
class Elem;
class Point;
class Node;
class EquationSystems;
class TiberLinearSolver;


template<typename T> class DenseMatrix;
template<typename T> class NumericVector;
template<typename T> class SparseMatrix;


class TBDLLOCAL MasterEquations : public SimulationInterface
{

  public:

    enum Solutions
    {
      eDensity,        /*!< electron density */
      hDensity,        /*!< hole density */
      Ec_edge,         /*!< conduction band lowest edge */
      Ev_edge          /*!< valence band highest edge */
      //Potential,        /*!< the potential */
      //ElField,            /*!< the field (negative gradient of potential) */
      //Displacement,     /*!< the electric displacement */
      //Polarization,
      //ChargeDensity     /*!< the source (charge density) */
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

        friend class MasterEquations;

    };


    // Destructor
    virtual ~MasterEquations(void);


    //! Create an DriftDiffusion object
    static MasterEquations* create(const ModelOptions& options);



    //! We need to create a physical model
    virtual PhysicalModel*
      create_bulk_model(const ModelOptions& options,
        const Material* mat) const;


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



  protected:

    //! The constructor
    /*!
     * Being private disables further inheritance.
     */
    MasterEquations(const ModelOptions& options);


    //virtual void get_conduction_band_edge(const Elem* elem, const std::vector<Point> pt);


    //virtual void get_valence_band_edge(const Elem* elem, const std::vector<Point>& pt);


    //! The initialization
    virtual void do_init(void);


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
        const std::vector<Elem>& elements);


    virtual void get_solution_secure(std::map<ID, std::vector<double> >& solutions);




  private:


    //! A static pointer to this
    static MasterEquations* _this;

    Device* _device;


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

   static void assemble(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
          SparseMatrix<Number>* jacobian);


    static void assemble_system(const NumericVector<Number>& x,
        NumericVector<Number>* residual,
          SparseMatrix<Number>* jacobian);


    //! The real assembly function
    //void do_assembly(EquationSystems& es, const std::string& system_name);



    template <int T>
    void do_assembly(const NumericVector<Number>& x, NumericVector<Number>* residual,
        SparseMatrix<Number>* jacobian);


    //! Rebuild the equation system if needed
    void rebuild_equation_system(void);



    //TiberLinearSystem *_sys_EcEv;

};


inline
MasterEquations*
MasterEquations::create(const ModelOptions& options)
{
  return new MasterEquations(options);
}


inline
MasterEquations::Options&
MasterEquations::get_my_options(void)
{
  return _options;
}


inline
MasterEquations::Options&
MasterEquations::get_options(void)
{
  return _options;
}


inline
unsigned int
MasterEquations::get_n_nonlinear_iterations(void) const
{
  return _n_nonlinear_iterations;
}

inline
double
MasterEquations::get_final_residual(void) const
{
  return _final_residual;
}



#endif // _MASTEREQUATIONS_H_
