// $Id: DefaultBandProperties.h 4184 2015-12-07 12:28:44Z maufder $

#ifndef _DEFAULTBANDPROPERTIES_H_
#define _DEFAULTBANDPROPERTIES_H_

#include "CarrierProperties.h"


//! Base class for band parameter models
class DefaultBandProperties : public CarrierProperties
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
