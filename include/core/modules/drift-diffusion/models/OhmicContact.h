// $Id$

#ifndef _OHMICCONTACT_H_
#define _OHMICCONTACT_H_

#include "ElectricalContact.h"


class TBDLLOCAL OhmicContact : public ElectricalContact
{
  public:

    //! Create an ohmic contact
    static OhmicContact* create(const ModelOptions& options);

  protected:

    //! The constructor
    OhmicContact(const ModelOptions& options);

    //! Create a new object
    virtual PhysicalModel* create_new(void) const;


    //! Calculate all coefficients
    virtual void do_compute(void);


};


//
// inline
//

inline
OhmicContact*
OhmicContact::create(const ModelOptions& options)
{
  return new OhmicContact(options);
}


inline
PhysicalModel*
OhmicContact::create_new(void) const
{
  return new OhmicContact(get_options());
}




#endif // _OHMICCONTACT_H_
