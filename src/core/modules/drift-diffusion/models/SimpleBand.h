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

    //! Read database
    virtual void read_database(void);

    //! Initialize
    virtual void do_init(void);


  private:

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
