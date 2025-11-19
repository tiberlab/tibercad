// $Id: DefaultBandProperties.C 4184 2015-12-07 12:28:44Z maufder $

#include "DefaultBandProperties.h"
#include "tibercad/base/ModelOptions.h"

#include "tibercad/module/TiberModule.h"


using namespace std;



DefaultBandProperties::DefaultBandProperties(const ModelOptions& options) :
    CarrierProperties(options)
{
}


DefaultBandProperties::~DefaultBandProperties(void)
{
}



DefaultBandProperties*
DefaultBandProperties::create(const ModelOptions& options)
{
  return new DefaultBandProperties(options);
}



