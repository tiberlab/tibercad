#include "DriftDiffusionHeatSource.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"




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


     _source_legend.resize(5);
     _source_legend[0]="Ejoule";
     _source_legend[1]="Hjoule";
     _source_legend[2]="RecSRH";
     _source_legend[3]="EPelTh";
     _source_legend[4]="HPelTh";

     _flux_legend.resize(2);
     _flux_legend[0]="Wn";
     _flux_legend[1]="Wp";
 
     

}




void
DriftDiffusionHeatSource::get_heat_sources(std::vector<Point> h_point,const Elem* elem,
					  std::vector< std::vector<double > >& heat_source)          
{
  // std::cout<<"start"<<std::endl;
  heat_source.clear();
  heat_source.resize(h_point.size());
  for(unsigned int n =0 ; n<h_point.size();n++)
  {
    heat_source[n].clear();
    heat_source[n].resize(5,0.0);
  }

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
      
      heat_source[n].resize(5);
      heat_source[n][0] = Ejoule;
      heat_source[n][1] = Hjoule;   
      heat_source[n][2] = Constants::e * R * (phih-phie + T * (Pp - Pn) ) ;
      heat_source[n][3] = pte;  
      heat_source[n][4] = pth;      
      
      
      
    } 
  }
  
   
    
    // for (unsigned int i=0; i<heat_source.size(); ++i)
    // std::cout<<heat_source[i][1]<<std::endl;

    //   std::cout<<h_point.size()<<std::endl;

  
}





void
DriftDiffusionHeatSource::get_power_fluxes(std::vector<Point> h_point,const Elem* elem,
					   std::vector<std::vector<RealGradient> >& power_fluxes,bool check_boundary)          
{
  
  power_fluxes.clear();
  power_fluxes.resize(h_point.size());
  for(unsigned int n =0 ; n<h_point.size();n++)
  {
    power_fluxes[n].clear();
    power_fluxes[n].resize(2);
  }
  
  
  
  std::vector< std::map< ID, double > > solution;
  
  if  (_simul->get_solution(elem,h_point,ID_set,solution) )
  {
    
    std::vector< std::map< ID, double > > solution;
    
    _simul->get_solution(elem,h_point,ID_set,solution); 
    
    for(unsigned int n =0 ; n<h_point.size();n++)
    {
      
      
      double Wn_x = solution[n].find(ID_vector[WNX])->second;
      double Wn_y = solution[n].find(ID_vector[WNY])->second;
      double Wn_z = solution[n].find(ID_vector[WNZ])->second; 
      
      
      double Wp_x = solution[n].find(ID_vector[WPX])->second;
      double Wp_y = solution[n].find(ID_vector[WPY])->second;
      double Wp_z = solution[n].find(ID_vector[WPZ])->second;
      
      power_fluxes[n][0](0) = Wn_x; 
      power_fluxes[n][0](1) = Wn_y; 
      power_fluxes[n][0](2) = Wn_z; 
      
      power_fluxes[n][1](0) = Wp_x; 
      power_fluxes[n][1](1) = Wp_y; 
      power_fluxes[n][1](2) = Wp_z;

    
    }

    
  }

 
}


