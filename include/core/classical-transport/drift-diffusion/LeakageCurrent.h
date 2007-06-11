// $Id$

#ifndef _LEAKAGECURRENT_H_
#define _LEAKAGECURRENT_H_


#include "ElectricalContact.h"
#include "Constants.h"
#include "DriftDiffusionProperties.h"


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


  protected:

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


  private:

    double _c;
    double _A;
    std::string _outer_voltage;


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
    _A(1e-21)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
}


inline
void
LeakageCurrent::do_init(void)
{
  ElectricalContact::do_init();

  _A = get_options().get_option("A", _A);
  _c = get_options().get_option("c", _c);
  _outer_voltage = get_options().get_option("outer_voltage", _outer_voltage);
}


inline
void
LeakageCurrent::get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c)
{
  a = 0.0;
  c = 0.0;

  double Vg = Variable::get_variable_value(_outer_voltage);
  double Ef = get_material().get_electron_electro_chemical_potential();
  double Vdiff = std::abs(Vg - Ef);

  double I = _A * Vdiff * Vdiff * std::exp(Vdiff / _c);
  if (Vg > Ef)
    I = -I;

  if (variable == DriftDiffusionDefs::FERMIE)
    c = - I / Constants::e;
}



#endif // _LEAKAGECURRENT_H_
