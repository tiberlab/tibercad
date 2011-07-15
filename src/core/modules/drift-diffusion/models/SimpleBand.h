// $Id$

#ifndef _SIMPLEBAND_H_
#define _SIMPLEBAND_H_


#include "BandProperties.h"


//! An implementation of BandProperties with explicit parameter values
class SimpleBand : public BandProperties
{

  public:

    //! Create a SimpleBand object
    static SimpleBand* create(const ModelOptions& options);


  protected:

    //! Constructor
    SimpleBand(const ModelOptions& options);

    //! Destructor
    virtual ~SimpleBand(void);

    //! We need some alloy parameters
    virtual void read_database_alloy(void);

    //! Read database
    virtual void read_database(void);

    //! Initialize
    virtual void do_init(void);

    //! Calculate band edges
    virtual void do_calculate(void);


  private:

    //! Reference energy (valence band edge)
    double _reference_energy;

    //! Band gap at 0 K
    double _bandgap;

    //! The bowing for Eg for alloys
    double _bow_Eg;

    //! Varshni parameters for the band gap
    double _varshni[2];

    //! The effective DOS, if given
    double _eff_DOS;

    //! Adjust DOS mass with degeneracy factor
    void _set_mdos(double& mdos);

    //! Set DOS mass from Nc
    void _set_mdos_from_Nc(double& Nc);

};


inline
SimpleBand*
SimpleBand::create(const ModelOptions& options)
{
  return new SimpleBand(options);
}

#endif // _SIMPLEBAND_H_
