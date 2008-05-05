#include "DriftDiffusionHeatSource.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"
#include "SimulationEnvironment.h"
#include "HeatModel.h"
#include "BoundaryProperties.h"
#include "Boundary.h"
//-------------------------------------------------------------------------//


//class  DriftDiffusionHeatSource::HeatSourceParameters;


 //! This class handles the heat source
 DriftDiffusionHeatSource::HeatSourceParameters::HeatSourceParameters(void):
    model_name("DriftDiffusion dissipation"),
    dd_simul_name("dd")
  {
    // TODO read default values from some text file
  }


void  DriftDiffusionHeatSource::copy_from(const PhysicalModelInterface *rhs)
{
  const DriftDiffusionHeatSource* mod = dynamic_cast<const DriftDiffusionHeatSource*> (rhs);
}

//-------------------------------------------------------------------------//


void   DriftDiffusionHeatSource::calculate_VCA (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const DriftDiffusionHeatSource* modA = dynamic_cast<const  DriftDiffusionHeatSource*>(comp_A);

  const DriftDiffusionHeatSource* modB = dynamic_cast<const  DriftDiffusionHeatSource*>(comp_B);


  //alloy(_kappa_e,modA->_kappa_e, modB->_kappa_e, xa);  

  // alloy(_kappa_h,modA->_kappa_h, modB->_kappa_h, xa); 
  
}


//---------------------------------------------------------//

void  DriftDiffusionHeatSource::do_init(void)
{


  heat_source_opt.model_name = "Drift diffusion dissipation";
  
  heat_source_opt.dd_simul_name = get_options().get_option("drift_diffusion_simulation", "no_current");

  _simul = SimulationInterface::find_simulation(heat_source_opt.dd_simul_name);

  
  if ( _simul == NULL)
   throw InitFailedException("Could not find " + heat_source_opt.dd_simul_name);


  std::string name_model = get_options(). get_option("recombination_type", "srh");
  std::string rec_string = "recHeat." + name_model;
 
 
  var_map.clear();
  var_map[EJOULE]=_simul->get_variable_id("HJouleN");

  var_map[HJOULE]=_simul->get_variable_id("HJouleP"); 
  //var_map[RECHEAT]=_simul->get_variable_id(rec_string); 
  var_map[RECHEAT]=_simul->get_variable_id("HRecomb");
  var_map[EPELTH]=_simul->get_variable_id("HPelThomE");  
  var_map[HPELTH]=_simul->get_variable_id("HPelThomH"); 
  var_map[WNX]=_simul->get_variable_id("PowerNx"); 
  var_map[WNY]=_simul->get_variable_id("PowerNy"); 
  var_map[WNZ]=_simul->get_variable_id("PowerNz");  
  var_map[WPX]=_simul->get_variable_id("PowerPx"); 
  var_map[WPY]=_simul->get_variable_id("PowerPy"); 
  var_map[WPZ]=_simul->get_variable_id("PowerPz");  


   std::map<ID,ID>::iterator      it(var_map.begin());
   std::map<ID,ID>::iterator      end(var_map.end());
   for(; it!=end; ++it)
     ID_set.insert(it->second);

}




void
DriftDiffusionHeatSource::get_heat_sources(std::vector<Point> h_point, const std::set<ID>& ids,
					   std::vector<std::map<ID, double> >& heat_sources)          

{
 
  
  heat_sources.clear();
  heat_sources.resize(h_point.size());
  for(unsigned int n =0 ; n<h_point.size();n++)
  {
    heat_sources[n].clear();
    //  heat_sources[n].resize(ids.size(),0.0);
  }
  
  const Elem*  elem = _heat_model->get_element();

  std::vector< std::map< ID, double > > solution;
 

  if  (_simul->get_solution(elem,h_point,ID_set,solution))
  {
   
    for(unsigned n =0; n<h_point.size();++n)
    { 
      
      

      double eJoule = solution[n].find(var_map[EJOULE])->second; 
      double hJoule = solution[n].find(var_map[HJOULE])->second;
      double RecHeat = solution[n].find(var_map[RECHEAT])->second;
      double ePelTh = solution[n].find(var_map[EPELTH])->second;
      double hPelTh = solution[n].find(var_map[HPELTH])->second;


      if  (ids.count(EJOULE))
	heat_sources[n][EJOULE]=eJoule;

      if  (ids.count(HJOULE))
	heat_sources[n][HJOULE]=hJoule;

      if  (ids.count(RECHEAT))
	heat_sources[n][RECHEAT]=RecHeat;

      if  (ids.count(EPELTH))
	heat_sources[n][EPELTH]=ePelTh;

      if  (ids.count(HPELTH))
	heat_sources[n][HPELTH]=hPelTh;

      if  (ids.count(100))
	heat_sources[n][100]=eJoule + hJoule + RecHeat + ePelTh + hPelTh;
        

    } 
  }
  
   
  
}


void
DriftDiffusionHeatSource::get_power_fluxes(std::vector<Point> h_point, const std::set<ID>& ids,
                                           std::vector<std::map<ID,RealGradient> >& power_fluxes)          
{
  
  power_fluxes.clear();
  power_fluxes.resize(h_point.size());
  for(unsigned int n =0 ; n<h_point.size();n++)
    power_fluxes[n].clear();
 

 
  const Elem*  elem = _heat_model->get_element();

  int side = _heat_model->get_side();
   //if side = -1 the check has not be done


  //if  no_check = true the check boundary is off
  bool do_calc = true;
  
  if (side >= 0 )
  { 
   
    const ElementSide elside(elem->top_parent(), side);
    
    do_calc =false; 
    
    if (_simul->get_environment().is_outer_boundary(elside))
      do_calc = true;

  }
      
  std::vector< std::map< ID, double > > solution;

  if  (_simul->get_solution(elem,h_point,ID_set,solution) && do_calc)
  {  
 
    for(unsigned int n =0 ; n<h_point.size();n++)
    {
     
      
          double Wn_x = solution[n].find(var_map[WNX])->second;
	  double Wn_y = solution[n].find(var_map[WNY])->second;
	  double Wn_z = solution[n].find(var_map[WNZ])->second; 
	  
	  double Wp_x = solution[n].find(var_map[WPX])->second;
	  double Wp_y = solution[n].find(var_map[WPY])->second;
	  double Wp_z = solution[n].find(var_map[WPZ])->second;
  
	  if (ids.count(0))
	  {
	    power_fluxes[n][0](0) = Wn_x; 
	    power_fluxes[n][0](1) = Wn_y; 
	    power_fluxes[n][0](2) = Wn_z; 
	  }

	  if (ids.count(1))
	  {
	    power_fluxes[n][1](0) = Wp_x; 
	    power_fluxes[n][1](1) = Wp_y; 
	    power_fluxes[n][1](2) = Wp_z;  
	  }
	  
	  if (ids.count(100))
	  {
	   
	    power_fluxes[n][100](0) = Wn_x + Wp_x; 
	    power_fluxes[n][100](1) = Wn_y + Wp_y; 
	    power_fluxes[n][100](2) = Wn_z + Wp_z;  
	   
	  }           
	  

    }

    
  }

}


std::map<ID,std::string>
DriftDiffusionHeatSource::get_source_legend(const std::set<std::string>& variables)
{


  if (variables.count("eJoule")     ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[EJOULE]="eJoule";

  
  
  if (variables.count("hJoule")     ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[HJOULE]="hJoule";
 
  
  if (variables.count("RecHeat")    ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[RECHEAT]="RecHeat";
  
  
  if (variables.count("ePelTh")     ||
      variables.count("hPelTh")     ||
      variables.count("thermal") )
    _source_legend[EPELTH]="ePelTh";
  
  if (variables.count("hPelTh")     ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[HPELTH]="hPelTh"; 
 

  return _source_legend;

}  


std::map<ID,std::string>
DriftDiffusionHeatSource::get_flux_legend(const std::set<std::string>& variables)
{

  if (variables.count("Wn")        ||
      variables.count("PowerFlux") ||
      variables.count("thermal") )
    _flux_legend[0]="Wn";

  if (variables.count("Wp")     ||
      variables.count("PowerFlux") ||
      variables.count("thermal") )
    _flux_legend[1]="Wp";

  
  return _flux_legend;

}  

void
DriftDiffusionHeatSource::do_print_info(void)
{
 
  std::string space = "           ";
  std::cout<<space<<"model:  "<<heat_source_opt.model_name<<std::endl;
  if (SimulationOptions::verbose() > 1)
    std::cout<<space<<"    DriftDiffusion simulation name:   "<<heat_source_opt.dd_simul_name<<std::endl;
  
}
