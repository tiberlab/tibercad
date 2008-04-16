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

     set_name("DD");

     std::string drift_diffusion_simulation = get_options().get_option("drift_diffusion_simulation", "no_current");

      _simul = SimulationInterface::find_simulation(drift_diffusion_simulation);

      if ( _simul == NULL)
	throw InitFailedException("Could not find " + drift_diffusion_simulation);

     std::string name_model = get_options(). get_option("recombination_type", "srh");
     std::string rec_string = "recombination." + name_model;




     ID_vector.resize(18);
     ID_vector[JNGRADPHIE]=_simul->get_variable_id("HJnGradPhie"); 
     ID_vector[JPGRADPHIH]=_simul->get_variable_id("HJpGradPhih"); 
     ID_vector[EJOULE]=_simul->get_variable_id("HJouleN");  
     ID_vector[HJOULE]=_simul->get_variable_id("HJouleP"); 
     ID_vector[PHIE]=_simul->get_variable_id("QFermi_e"); 
     ID_vector[PHIH]=_simul->get_variable_id("QFermi_h"); 
     ID_vector[PN]=_simul->get_variable_id("Pn"); 
     ID_vector[PP]=_simul->get_variable_id("Pp");  
     ID_vector[TEMP]=_simul->get_variable_id("Temp"); 
     ID_vector[SRHREC]=_simul->get_variable_id(rec_string); 
     ID_vector[WNX]=_simul->get_variable_id("PowerNx"); 
     ID_vector[WNY]=_simul->get_variable_id("PowerNy"); 
     ID_vector[WNZ]=_simul->get_variable_id("PowerNz");  
     ID_vector[WPX]=_simul->get_variable_id("PowerPx"); 
     ID_vector[WPY]=_simul->get_variable_id("PowerPy"); 
     ID_vector[WPZ]=_simul->get_variable_id("PowerPz");  
     ID_vector[PELTHE]=_simul->get_variable_id("HPelThomE"); 
     ID_vector[PELTHH]=_simul->get_variable_id("HPelThomH");  


     for (unsigned int n=0; n<ID_vector.size(); ++n)
       ID_set.insert(ID_vector[n]);


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
      
      double R      = solution[n].find(ID_vector[SRHREC]  )->second;
      double Ejoule = solution[n].find(ID_vector[EJOULE]  )->second; 
      double Hjoule = solution[n].find(ID_vector[HJOULE]  )->second; 
      double Pn     = solution[n].find(ID_vector[PN]      )->second; 
      double Pp     = solution[n].find(ID_vector[PP]      )->second; 
      double T      = solution[n].find(ID_vector[TEMP]    )->second; 
      double phie   = solution[n].find(ID_vector[PHIE]    )->second; 
      double phih   = solution[n].find(ID_vector[PHIH]    )->second; 
      double pte    = solution[n].find(ID_vector[PELTHE]  )->second; 
      double pth    = solution[n].find(ID_vector[PELTHH]  )->second; 
      
      double RecHeat = Constants::e * R * (phih-phie + T * (Pp - Pn) );
      

      if  (ids.count(0))
	heat_sources[n][0]=Ejoule;
      if  (ids.count(1))
	heat_sources[n][1]=Hjoule;
      if  (ids.count(2))
	heat_sources[n][2]=RecHeat;
      if  (ids.count(3))
	heat_sources[n][3]=pte;
      if  (ids.count(4))
	heat_sources[n][4]=pth;
      if  (ids.count(100))
	heat_sources[n][5]=Ejoule + Hjoule + RecHeat + pte + pth;
        

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
  {
    power_fluxes[n].clear();
    //    power_fluxes[n].resize(ids.size()); 
  }
  
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
     
	  double Wn_x = solution[n].find(ID_vector[WNX])->second;
	  double Wn_y = solution[n].find(ID_vector[WNY])->second;
	  double Wn_z = solution[n].find(ID_vector[WNZ])->second; 
	  
	  double Wp_x = solution[n].find(ID_vector[WPX])->second;
	  double Wp_y = solution[n].find(ID_vector[WPY])->second;
	  double Wp_z = solution[n].find(ID_vector[WPZ])->second;
      

  
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
	    power_fluxes[n][1](0) = Wn_x + Wp_x; 
	    power_fluxes[n][1](1) = Wn_y + Wp_y; 
	    power_fluxes[n][1](2) = Wn_z + Wp_z;  
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
    _source_legend[0]="eJoule";

  
  
  if (variables.count("hJoule")     ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[1]="hJoule";
 
  
  if (variables.count("RecHeat")    ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[2]="RecHeat";
  
  
  if (variables.count("ePelTh")     ||
      variables.count("hPelTh")     ||
      variables.count("thermal") )
    _source_legend[3]="ePelTh";
  
  if (variables.count("hPelTh")     ||
      variables.count("HeatSource") ||
      variables.count("thermal") )
    _source_legend[4]="hPelTh"; 
 

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
