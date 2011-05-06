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


    //! Calculate all coefficients
    virtual void do_compute(void);


  private:

    //! The metal Fermi level (negative of work function)
    double _metal_Ef;

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






#endif // _SCHOTTKYCONTACT_H_
