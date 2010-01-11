// $Id$

#ifndef _SCHOTTKYCONTACT_H_
#define _SCHOTTKYCONTACT_H_

#include "ElectricalContact.h"



class SchottkyContact : public ElectricalContact
{
  public:

    //! Create a schottky contact
    static SchottkyContact* create(const ModelOptions& options);

    //! \copydoc ElectricalContact::get_boundary_value()
    virtual double get_boundary_value(DriftDiffusionDefs::DDVariable variable);


  protected:

    //! The constructor
    SchottkyContact(const ModelOptions& options);

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
SchottkyContact::create(const ModelOptions& options)
{
  return new SchottkyContact(options);
}



inline
SchottkyContact::SchottkyContact(const ModelOptions& options)
  : ElectricalContact(options),
    _band('c'),
    _fixed_barrier(true)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::DIRICHLET);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::DIRICHLET);
}




#endif // _SCHOTTKYCONTACT_H_
