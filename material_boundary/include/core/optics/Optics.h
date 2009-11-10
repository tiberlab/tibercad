// $Id$

#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "SimulationInterface.h"


//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point
 */
class Optics : public SimulationInterface
{

  public:

    //! The constructor
    Optics(void) { };

    //! The destructor
    virtual ~Optics(void) { };


  protected:

    //! do_solve() needs to be implemented by derived classes
    virtual void do_solve(void) = 0;


  private:

};


#endif // _OPTICS_H_
