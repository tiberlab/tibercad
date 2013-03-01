// $Id: ConstantHeatSource.C 2457 2011-03-06 23:52:12Z gromano $

#include "ConstantHeatSource.h"
#include "Material.h"

#include "TiberModule.h"




using namespace std;


ConstantHeatSource::ConstantHeatSource(const ModelOptions& options):HeatSourceModel(options)
{
}

void
ConstantHeatSource::do_init(void)
{
  double heat_source = 0.0;
 
  get_parameter("H",heat_source);  

  set_heat_source(heat_source);

}




 
