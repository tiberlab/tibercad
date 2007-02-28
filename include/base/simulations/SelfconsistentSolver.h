// $Id$

#ifndef _SELFCONSISTENTSOLVER_H_
#define _SELFCONSISTENTSOLVER_H_

#include "SimulationInterface.h"


//! Interface for self-consistent calculations
class SelfconsistentSolver : public SimulationInterface
{

  public:

    //! The destructor
    virtual ~SelfconsistentSolver(void);

    //! Create a Sweep object
    static SelfconsistentSolver* create(void);


  protected:

    //! The empty Constructor
    SelfconsistentSolver(void);

    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);

    /*! \copydoc SimulationInterface::do_equilibrium() */
    virtual void do_equilibrium(void);

    /*! \copydoc SimulationInterface::do_solve() */
    virtual void do_solve(void);

    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);
    

  private:

    //! The first simulation
    SimulationInterface* _simulation1;
    
    //! The second simulation
    SimulationInterface* _simulation2;

    //! The maximum number of iterations
    unsigned int _max_it;

    //! The relative tolerance
    double _rel_tol;

    //! The absolute tolerance
    double _abs_tol;

    //! The relaxation factor to be used
    double _relax;


    //! Get the norm of the difference of two solutions
    /*!
     * The two vectors have to be of the same size
     */
    double get_norm_of_difference(NumericVector<double>& vec1,
        NumericVector<double>& vec2);
    
};


//
// inline methods
//

inline
SelfconsistentSolver::SelfconsistentSolver(void)
  : _simulation1(NULL),
    _simulation2(NULL),
    _max_it(5),
    _rel_tol(1e-3),
    _abs_tol(1e-3),
    _relax(0.5)
{
}

inline
SelfconsistentSolver::~SelfconsistentSolver(void)
{
}

inline
SelfconsistentSolver*
SelfconsistentSolver::create(void)
{
  return new SelfconsistentSolver();
}



#endif // _SELFCONSISTENTSOLVER_H_
