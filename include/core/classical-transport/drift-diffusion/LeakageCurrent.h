// $Id$

#ifndef _LEAKAGECURRENT_H_
#define _LEAKAGECURRENT_H_


#include "ElectricalContact.h"
#include "Constants.h"


//! A simple class to impose a current density on a boundary
/*!
 * For now one can only impose electron current density, and in the
 * input file it is called voltage
 */
class LeakageCurrent : public ElectricalContact
{

  public:

    //! The constructor
    LeakageCurrent(void);
    
    //! The destructor
    ~LeakageCurrent(void) {};

    //! Create an object
    static LeakageCurrent* create(void);

    //! Get normal derivative
    virtual void get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c);

  private:


};


//
// inline members
// 

inline
LeakageCurrent*
LeakageCurrent::create(void)
{
  return new LeakageCurrent();
}



inline
LeakageCurrent::LeakageCurrent(void)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
}


inline
void
LeakageCurrent::get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c)
{
  a = 0.0;
  c = 0.0;

  if (variable == DriftDiffusionDefs::FERMIE)
    c = - get_variable_value() / Constants::e;
}



#endif // _LEAKAGECURRENT_H_
