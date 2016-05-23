
#ifndef _GUESTINTERFACE_H_
#define _GUESTINTERFACE_H_


#include "ElectricalContact.h"


class TBDLLOCAL GuestInterface : public ElectricalContact
{
  public:

    //! Create a schottky contact
    static GuestInterface* create(const ModelOptions& options);


  protected:

    //! The contructor
    GuestInterface(const ModelOptions& options);

    //! Initialize
    virtual void do_init(void);


    //! Setup the coefficients
    virtual void do_compute(void);


    //! Get the contact voltage drop
    /*!
     * A positive value means: the inner voltage is lower than
     * the outer contact voltage
     */
    //double get_contact_voltage_drop() const;


  private:

  ID _host_sim;
};

inline
GuestInterface*
GuestInterface::create(const ModelOptions& options)
{
  return new GuestInterface(options);
}



#endif // _GUESTINTERFACE_H_