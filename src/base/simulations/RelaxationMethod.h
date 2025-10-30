// $Id$

#ifndef _RELAXATIONMETHOD_H_
#define _RELAXATIONMETHOD_H_

#include "SelfconsistentSolver.h"


//! Interface for self-consistent calculations
class TBDLLOCAL RelaxationMethod : public SelfconsistentSolver
{

  public:

    //! The destructor
    virtual ~RelaxationMethod(void);

    //! Create a Sweep object
    static RelaxationMethod* create(const ModelOptions& options);


  protected:

    //! The empty Constructor
    RelaxationMethod(const ModelOptions& options);

    
    /*! \copydoc SimulationInterface::do_solve() */
    virtual void do_solve(void);

    
    /*! \copydoc SimulationInterface::parse_options() */
    virtual void parse_options(void);



  private:

    //! The relaxation factor to be used
    double _relax;


};


//
// inline methods
//

inline
RelaxationMethod::RelaxationMethod(const ModelOptions& options)
  : SelfconsistentSolver(options),
    _relax(1.0)
{
}

inline
RelaxationMethod::~RelaxationMethod(void)
{
}

inline
RelaxationMethod*
RelaxationMethod::create(const ModelOptions& options)
{
  return new RelaxationMethod(options);
}



#endif // _RELAXATIONMETHOD_H_
