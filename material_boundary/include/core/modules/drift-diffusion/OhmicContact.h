// $Id$

#ifndef _OHMICCONTACT_H_
#define _OHMICCONTACT_H_

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"


class OhmicContact : public ElectricalContact
{
  public:

    //! The constructor
    OhmicContact(void);

    //! Create an ohmic contact
    static OhmicContact* create(void);

    //! \copydoc ElectricalContact::get_boundary_value()
    virtual double get_boundary_value(DriftDiffusionDefs::DDVariable variable);
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
OhmicContact::get_boundary_value(DriftDiffusionDefs::DDVariable variable)
{
  double val = 0.0;
  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      val = get_reference_material().get_equilibrium_fermi_level();
      break;
    case DriftDiffusionDefs::FERMIE:
      break;
    case DriftDiffusionDefs::FERMIH:
      break;
    /*
    case DriftDiffusionDefs::DENSE:
      {
        DriftDiffusionProperties& sc = get_material();
        sc.set_potentials(get_material().get_equilibrium_fermi_level());
        sc.calculate_densities();
        val = sc.get_electron_density();
      }
      break;
    case DriftDiffusionDefs::DENSH:
      {
        DriftDiffusionProperties& sc = get_material();
        sc.set_potentials(get_material().get_equilibrium_fermi_level());
        sc.calculate_densities();
        val = sc.get_hole_density();
      }
      break;
    */
  }

  return val;
}

#endif // _OHMICCONTACT_H_
