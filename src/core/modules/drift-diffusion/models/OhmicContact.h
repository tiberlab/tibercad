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






#endif // _OHMICCONTACT_H_
