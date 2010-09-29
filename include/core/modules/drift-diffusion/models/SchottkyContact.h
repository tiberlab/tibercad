// $Id$

#ifndef _SCHOTTKYCONTACT_H_
#define _SCHOTTKYCONTACT_H_

#include "ElectricalContact.h"


/*!
 * \brief A Schottky contact
 */
class TBDLLOCAL SchottkyContact : public ElectricalContact
{
  public:

    //! Create a schottky contact
    static SchottkyContact* create(const ModelOptions& options);


  protected:

    //! The constructor
    SchottkyContact(const ModelOptions& options);

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);

    //! Create a new one
    virtual PhysicalModel* create_new(void) const;


    //! Calculate all coefficients
    virtual void do_compute(void);


  private:

    //! The work function
    double _workfunction;

    //! The reference band
    char _band;

    //! Is this a fixed barrier or not?
    bool _fixed_barrier;


    //! Do we include thermionic emission?
    bool _thermionic_emission;

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
PhysicalModel*
SchottkyContact::create_new(void) const
{
  return new SchottkyContact(get_options());
}





#endif // _SCHOTTKYCONTACT_H_
