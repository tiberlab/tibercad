// $Id$

#ifndef _OHMICCONTACT_H_
#define _OHMICCONTACT_H_

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"

using namespace DriftDiffusionDefs;

class OhmicContact : public ElectricalContact
{
  public:

    //! The constructor
    OhmicContact(void);

    //! Create an ohmic contact
    static OhmicContact* create(void);
    
    //! \copydoc ElectricalContact::get_boundary_value()
    virtual double get_boundary_value(DriftDiffusionDefs::Variable variable);
};


//
// inline
// 

inline
OhmicContact*
OhmicContact::create(void)
{
  return new OhmicContact();
}


inline
OhmicContact::OhmicContact(void)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::DIRICHLET);
}

inline
double
OhmicContact::get_boundary_value(DriftDiffusionDefs::Variable variable)
{
  double val = 0.0;
  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      val = get_material().get_equilibrium_fermi_level();
      break;
    case DriftDiffusionDefs::FERMIE:
      break;
    case DriftDiffusionDefs::FERMIH:
      break;
  }
  
  return val;
}

#endif // _OHMICCONTACT_H_
