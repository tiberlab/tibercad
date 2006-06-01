#ifndef _SCHOTTKYCONTACT_H_
#define _SCHOTTKYCONTACT_H_

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"

using namespace DriftDiffusionDefs;

class SchottkyContact : public ElectricalContact
{
  public:

    //! The constructor
    SchottkyContact(const std::string identifier);
    
    //! \copydoc ElectricalContact::get_boundary_value()
    virtual double get_boundary_value(DriftDiffusionDefs::Variable variable);

    //! Set the Schottky barrier
    void set_schottky_barrier(double barrier);

  private:

    //! The Schottky barrier
    double _barrier;
};

inline
SchottkyContact::SchottkyContact(const std::string identifier)
  : ElectricalContact(identifier)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::DIRICHLET);
}

inline
double
SchottkyContact::get_boundary_value(DriftDiffusionDefs::Variable variable)
{
  double val = 0.0;
  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      val = get_material().get_conduction_band_edge() - _barrier;
      break;
    case DriftDiffusionDefs::FERMIE:
      break;
    case DriftDiffusionDefs::FERMIH:
      break;
  }
  
  return val;
}

inline
void
SchottkyContact::set_schottky_barrier(double barrier)
{
  _barrier = barrier;
}

#endif // _SCHOTTKYCONTACT_H_
