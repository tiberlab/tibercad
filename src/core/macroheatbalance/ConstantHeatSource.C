// $Id$

#include "ConstantHeatSource.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"
#include "SimulationEnvironment.h"
#include "SimulationOptions.h"
#include "HeatModel.h"
#include "BoundaryProperties.h"
#include "Boundary.h"
//-------------------------------------------------------------------------//



//-------------------------------------------------------------------------//


void  ConstantHeatSource::do_init_alloy (const PhysicalModelInterface *comp_A,
                                                const PhysicalModelInterface *comp_B, double xa)
{
  const ConstantHeatSource* modA = dynamic_cast<const  ConstantHeatSource*>(comp_A);

  const ConstantHeatSource* modB = dynamic_cast<const  ConstantHeatSource*>(comp_B);

}


//---------------------------------------------------------//

void  ConstantHeatSource::do_init(void)
{

  _heat_source = get_options().get_option("H",0.0);
 

}





void
ConstantHeatSource::get_power_fluxes(const Elem*  elem,std::vector<Point> h_point,std::vector<std::map<ID,RealGradient> >& power_fluxes)
{

  power_fluxes.clear();
  power_fluxes.resize(h_point.size());

  for(unsigned int n =0 ; n<h_point.size();n++)
    power_fluxes[n].clear();

}
void
ConstantHeatSource::get_heat_sources(const Elem*  elem,std::vector<Point> h_point, std::vector<std::map<ID, double> >& heat_sources)

{


  heat_sources.clear();
  heat_sources.resize(h_point.size());
 

  for (ID n =0;n<h_point.size();n++)
    heat_sources[n][0]=_heat_source;
   
     
  
}

std::map<ID,std::string> ConstantHeatSource::get_source_legend(const std::set<std::string>& variables)
{


  if (variables.count("HeatSource"))
      _source_legend[0]="H";

  return  _source_legend;

}

void
ConstantHeatSource::do_print_info(void)
{  
  if (SimulationOptions::verbose() > 1)
  {
    std::string space = "           ";
    std::cout<<space<<"model:  "<<get_name()<<std::endl;
    std::cout<<space<<"    HeatSource:   "<<_heat_source<<" W/cm3  "<<std::endl;
  }

}
