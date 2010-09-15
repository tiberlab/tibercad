// $Id$

#include "ConstantHeatSource.h"
#include "Material.h"


// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(ConstantHeatSource, heat_source, constant)

using namespace std;


ConstantHeatSource::ConstantHeatSource(const ModelOptions& options):HeatSourceModel(options)
{
}

void
ConstantHeatSource::do_init(void)
{
  double heat_source = 0.0;
  get_parameter("H", heat_source,true);
  set_heat_source(heat_source);

}




 
