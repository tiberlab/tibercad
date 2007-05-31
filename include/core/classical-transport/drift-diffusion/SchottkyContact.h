// $Id$

#ifndef _SCHOTTKYCONTACT_H_
#define _SCHOTTKYCONTACT_H_

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"

using namespace DriftDiffusionDefs;

class SchottkyContact : public ElectricalContact
{
  public:

    //! The constructor
    SchottkyContact(void);

    //! Create a schottky contact
    static SchottkyContact* create(void);
    
    //! \copydoc ElectricalContact::get_boundary_value()
    virtual double get_boundary_value(DriftDiffusionDefs::Variable variable);


  protected:

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


  private:

    //! The Schottky barrier
    double _barrier;


    //! Whether to consider doping dependence of barrier height
    bool _doping_dependent;
};



//
// inline
// 


inline
SchottkyContact*
SchottkyContact::create(void)
{
  return new SchottkyContact();
}



inline
SchottkyContact::SchottkyContact(void)
  : _doping_dependent(false)
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
SchottkyContact::do_init(void)
{
  ElectricalContact::do_init();

  _doping_dependent = get_options().get_option("doping_dependent_height", false); 
  _barrier = get_options().get_option("barrier_height", 0.8);
}




#endif // _SCHOTTKYCONTACT_H_
