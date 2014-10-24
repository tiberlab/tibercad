// $Id$

#ifndef _KPBAND_H_
#define _KPBAND_H_


#include "BandProperties.h"

class DDsemiconductor;

//! An implementation of BandProperties with \f[k \cdot p\f] values
class KPBand : public BandProperties
{

  public:

    //! Create a KPBand object
    static KPBand* create(const ModelOptions& options);


  protected:

    //! Constructor
    KPBand(const ModelOptions& options);

    //! Destructor
    virtual ~KPBand(void);

    //! Create the kp bulk model
    virtual void prepare_submodels(void);

    //! Initialize
    virtual void do_init(void);

    //! Read database
    virtual void read_database(void);

    //! \copydoc PhysicalModelInterface::do_print_info(void)
    virtual void do_print_info(void);

    //! Calculate band properties
    virtual void do_calculate(void);


  private:

    //! The physical model for this semiconductor
    DDsemiconductor* _bulk_model;

    //! An alternative DOS mass, not from kp
    /*!
     * For holes, warping leads to errors in the DOS mass calculation
     */
    double _dos_mass;

};


inline
KPBand*
KPBand::create(const ModelOptions& options)
{
  return new KPBand(options);
}

#endif // _KPBAND_H_
