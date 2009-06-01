// $Id$

#ifndef _SCHOTTKYCONTACT_H_
#define _SCHOTTKYCONTACT_H_

#include "ElectricalContact.h"



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

    //! The work function
    double _workfunction;

    //! The reference band
    char _band;

    //! Is this a fixed barrier or not?
    bool _fixed_barrier;

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
  : _band('c'), _fixed_barrier(true)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::DIRICHLET);
}




#endif // _SCHOTTKYCONTACT_H_
