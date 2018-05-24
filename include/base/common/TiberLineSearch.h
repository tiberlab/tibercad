// $Id$


#ifndef _TIBERLINESEARCH_H_
#define _TIBERLINESEARCH_H_


#include "TiberNonlinearSystem.h"

class TiberLinearSolver;


//! A base class for the TiberCAD line search algorithms
class TBDLLOCAL TiberLineSearch : public TiberNonlinearSystem
{

  public:

    //! Destructor
    virtual ~TiberLineSearch(void);


    /*! \copydoc TiberNonlinearSystem::clear() */
    virtual void clear(void);


    /*! \copydoc TiberNonlinearSystem::reinit() */
    virtual void reinit(void);


    /*! \copydoc System:user_initialization() */
    virtual void user_initialization(void);


    //! Get the solution vector
    //virtual libMesh::NumericVector<double>& get_solution_vector(void);

    

  protected:

    //! Constructor
    TiberLineSearch(libMesh::EquationSystems& es,
        const std::string& name, const unsigned int number);


    //! Get the linear solver
    TiberLinearSolver* get_linear_solver(void);
 

    //! Get the relative nonlinear tolerance
    double get_nonlinear_rtol(void) const;

    //! Get the absolute nonlinear tolerance
    double get_nonlinear_atol(void) const;

    //! Get the nonlinear step tolerance
    double get_nonlinear_stol(void) const;

    //! Get the maximum number of nonlinear iterations
    int get_nonlinear_max_it(void) const;

    //! Get the line search max step (for l2 norm)
    double get_max_step(void) const;

    //! Get the line search max absolute step
    double get_max_abs_step(void) const;

    //! Get divergence tolerance
    double get_divergence_tol(void) const;

    /*! \copydoc TiberEqSystem::parse_options() */
    virtual void parse_options(void);



  private:

     //! The parent class type
    typedef TiberNonlinearSystem Parent;


    //! The linear solver to be used for the Newton iteration
    TiberLinearSolver* _solver;

    
    //! The relative linear tolerance
    double _nonlinear_rtol;

    //! The absolute linear tolerance
    double _nonlinear_atol;

    //! The step tolerance
    double _nonlinear_stol;
      
    //! The maximum number of nonlinear iterations
    int _nonlinear_max_it;

    //! The absolute maximum search step
    double _max_abs_step;

    //! The maximum search step
    double _max_step;

    //! The divergence tolerance
    /*!
     * If the step norm increases by a factor of \c _divergence_tol or more
     * we assume failure.
     */
    double _divergence_tol;



};


//
// inline methods
//


//inline
//libMesh::NumericVector<double>&
//TiberLineSearch::get_solution_vector(void)
//{
  //return get_vector("sol");
//  return *current_local_solution;
//}


inline
TiberLinearSolver*
TiberLineSearch::get_linear_solver(void)
{
  return _solver;
}


inline
double
TiberLineSearch::get_nonlinear_rtol(void) const
{
  return _nonlinear_rtol;
}


inline
double
TiberLineSearch::get_nonlinear_atol(void) const
{
  return _nonlinear_atol;
}


inline
double
TiberLineSearch::get_nonlinear_stol(void) const
{
  return _nonlinear_stol;
}



inline
int
TiberLineSearch::get_nonlinear_max_it(void) const
{
  return _nonlinear_max_it;
}


inline
double
TiberLineSearch::get_max_abs_step(void) const
{
  return _max_abs_step;
}




inline
double
TiberLineSearch::get_max_step(void) const
{
  return _max_step;
}



inline
double
TiberLineSearch::get_divergence_tol(void) const
{
  return _divergence_tol;
}



#endif // _TIBERLINESEARCH_H_
