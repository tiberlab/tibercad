// $Id$

#ifndef _LEAKAGECURRENT_H_
#define _LEAKAGECURRENT_H_


#include "ElectricalContact.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"

class TunnelingCurrent;

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
    virtual void get_normal_derivative(DriftDiffusionDefs::DDVariable variable,
        double& a, double& c);


  protected:

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


  private:

    double _c;
    double _A;
    std::string _outer_voltage;

    TunnelingCurrent* _tc;

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
  : _c(0.02585),
    _A(1e-21),
    _tc(NULL)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
}



#endif // _LEAKAGECURRENT_H_
