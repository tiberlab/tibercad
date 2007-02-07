// $Id$

#ifndef _FERMILEVELPINNING_H_
#define _FERMILEVELPINNING_H_

#include <iostream>
#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"

using namespace DriftDiffusionDefs;

class FermiLevelPinning : public ElectricalContact
{
  public:

    //! The constructor
    FermiLevelPinning(void);

    //! Create an ohmic contact
    static FermiLevelPinning* create(void);
    
    //! \copydoc ElectricalContact::get_boundary_value()
    virtual double get_boundary_value(DriftDiffusionDefs::Variable variable);

    //! set the value of pinning
    void set_pinning(double pinning);


  protected:

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


  private:

    //! The value of the pinning
    /*!
     * The pinning is defined by the potential
     * \f[E_c = E_c^0 - \varphi = \phi_n + \phi_pinning\f]
     */
    double _pinning;
};


//
// inline
// 

inline
FermiLevelPinning*
FermiLevelPinning::create(void)
{
  return new FermiLevelPinning();
}


inline
FermiLevelPinning::FermiLevelPinning(void)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::PINNING);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
}

inline
void
FermiLevelPinning::set_pinning(double pinning)
{
  _pinning = pinning;
}

inline
double
FermiLevelPinning::get_boundary_value(DriftDiffusionDefs::Variable variable)
{
  double val = 0.0;
  if (variable == DriftDiffusionDefs::POTENTIAL)
    val = get_material().get_conduction_band_edge() - _pinning;

  return val;
}


inline
void
FermiLevelPinning::do_init(void)
{
  ElectricalContact::do_init();

  set_pinning(get_options().get_option("pinning", 0.8));
}


#endif // _FERMILEVELPINNING_H_
