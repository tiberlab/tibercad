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

    //! The destructor
    ~LeakageCurrent(void) {};

    //! Create an object
    static LeakageCurrent* create(const ModelOptions& options);

    //! Get normal derivative
    virtual void get_normal_derivative(DriftDiffusionDefs::DDVariable variable,
        double& a, double& c);


  protected:

    //! The constructor
    LeakageCurrent(const ModelOptions& options);

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


  private:

    double _c_param;
    double _A_param;
    std::string _outer_voltage;

    TunnelingCurrent* _tc;

};


//
// inline members
//

inline
LeakageCurrent*
LeakageCurrent::create(const ModelOptions& options)
{
  return new LeakageCurrent(options);
}



inline
LeakageCurrent::LeakageCurrent(const ModelOptions& options)
  : ElectricalContact(options),
    _c_param(0.02585),
    _A_param(1e-21),
    _tc(NULL)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
}



#endif // _LEAKAGECURRENT_H_
