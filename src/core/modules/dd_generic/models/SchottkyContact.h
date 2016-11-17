// $Id: SchottkyContact.h 4135 2015-09-25 10:19:38Z maufder $

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
    std::string _band;

    //! Is this a fixed barrier or not?
    bool _fixed_barrier;


    //! Do we include thermionic emission?
    bool _thermionic_emission;

    //! Do we use Scott and Malliaras field and mobility dependent recombination velocity? (Chem. Phys. Lett. 299 (1999) 115)
    bool _scott_malliaras;

    //! Image-force lowering
    bool _barrier_lowering;
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
