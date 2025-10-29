// $Id$

#ifndef _DEFAULTBANDPROPERTIES_H_
#define _DEFAULTBANDPROPERTIES_H_

#include "BandProperties.h"


//! Base class for band parameter models
class DefaultBandProperties : public BandProperties
{

  public:

    //! Destructor
    virtual ~DefaultBandProperties(void); 

    //! Creator method
    static DefaultBandProperties* create(const ModelOptions& options);


  protected:

    //! Constructor
    DefaultBandProperties(const ModelOptions& options);

  private:

};


#endif // _DEFAULTBANDPROPERTIES_H_
