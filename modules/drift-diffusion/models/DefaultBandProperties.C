// $Id$

#include "DefaultBandProperties.h"
#include "tibercad/base/ModelOptions.h"

#include "tibercad/module/TiberModule.h"


using namespace std;



DefaultBandProperties::DefaultBandProperties(const ModelOptions& options) :
    BandProperties(options)
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



