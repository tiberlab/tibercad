#include "DriftDiffusionHeatSource.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"
#include "SimulationEnvironment.h"



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

     std::string name_model = get_options().get_option("recombination_type", "srh");
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
     ID_vector[TOLD]=_simul->get_variable_id("Temp"); 
     ID_vector[SRHREC]=_simul->get_variable_id(rec_string); 
     ID_vector[PELTIERNX]=_simul->get_variable_id("PeltierNx"); 
     ID_vector[PELTIERNY]=_simul->get_variable_id("PeltierNy"); 
     ID_vector[PELTIERNZ]=_simul->get_variable_id("PeltierNz");  
     ID_vector[PELTIERPX]=_simul->get_variable_id("PeltierPx"); 
     ID_vector[PELTIERPY]=_simul->get_variable_id("PeltierPy"); 
     ID_vector[PELTIERPZ]=_simul->get_variable_id("PeltierPz"); 
     ID_vector[PELTHOME]=_simul->get_variable_id("HPelThomE"); 
     ID_vector[PELTHOMH]=_simul->get_variable_id("HPelThomH"); 

     for (unsigned int n=0; n<ID_vector.size(); ++n)
       ID_set.insert(ID_vector[n]);
     

     _legend.resize(5);
     _legend[0]="Ejoule";
     _legend[1]="Hjoule";
     _legend[2]="Rec";
     _legend[3]="EPelTh";
     _legend[4]="HPelTh";


}



void
DriftDiffusionHeatSource::get_heat_source(std::vector<Point> h_point,const Elem* elem,
                                	   std::vector< double > & heat_source)          
{
  heat_source.clear();
  heat_source.resize(h_point.size());

  SimulationEnvironment& env =  _simul->get_environment();  
    
    if  (env.contains_element(elem) & _simul->is_initialized())
    {
      std::vector< std::map< ID, double > > solution;

      _simul->get_solution(elem,h_point,ID_set,solution); 

      for(unsigned n =0; n<h_point.size();++n)
      {

	double R     = solution[n].find(ID_vector[SRHREC])->second;
	double phi_h = solution[n].find(ID_vector[PHIH])->second;
	double phi_e = solution[n].find(ID_vector[PHIE])->second;

	double JnGradPhiE = solution[n].find(ID_vector[JNGRADPHIE])->second; 
	double JpGradPhiH = solution[n].find(ID_vector[JPGRADPHIH])->second;
        double SrhSource =  Constants::e * R * (phi_h - phi_e);

	heat_source[n]= JnGradPhiE + JpGradPhiH + SrhSource;

	
      }
      
    }

 
}

void
DriftDiffusionHeatSource::get_heat_source_output(std::vector<Point> h_point,const Elem* elem,
						 std::vector< std::vector<double > >& heat_source)          
{
   heat_source.clear();
   heat_source.resize(h_point.size());

  SimulationEnvironment& env =  _simul->get_environment();  
    
    if  (env.contains_element(elem) & _simul->is_initialized())
    {
      std::vector< std::map< ID, double > > solution;
 
      _simul->get_solution(elem,h_point,ID_set,solution); 

      for(unsigned n =0; n<h_point.size();++n)
      {

	double R      = solution[n].find(ID_vector[SRHREC])->second;
        double Ejoule = solution[n].find(ID_vector[EJOULE])->second; 
        double Hjoule = solution[n].find(ID_vector[HJOULE])->second; 
	double Pn     = solution[n].find(ID_vector[PN]    )->second; 
	double Pp     = solution[n].find(ID_vector[PP]    )->second; 
        double T      = solution[n].find(ID_vector[TOLD]  )->second; 
	double phie   = solution[n].find(ID_vector[PHIE]  )->second; 
	double phih   = solution[n].find(ID_vector[PHIH]  )->second; 
	double pte   = solution[n].find(ID_vector[PELTHOME]  )->second; 
	double pth   = solution[n].find(ID_vector[PELTHOMH]  )->second; 
       
        
	heat_source[n].push_back(Ejoule);
	heat_source[n].push_back(Hjoule);   
	heat_source[n].push_back(Constants::e * R * (phih-phie + T * (Pp - Pn) ) );
	// heat_source[n].push_back(Constants::e * R * T * (Pp - Pn) ) ;
       	heat_source[n].push_back(pte);  
	heat_source[n].push_back(pth);      

	
       } 
      
    }

}

void
DriftDiffusionHeatSource::get_flux_heat_source(std::vector<Point> h_point,const Elem* elem,
					       std::vector<RealGradient>& vector_heat_source)          
{

  vector_heat_source.clear();
  vector_heat_source.resize(h_point.size());

  SimulationEnvironment& env =  _simul->get_environment();  
    
  if  (env.contains_element(elem) & _simul->is_initialized())
  {

    std::vector< std::map< ID, double > > solution;

     _simul->get_solution(elem,h_point,ID_set,solution); 

    for(unsigned int n =0 ; n<h_point.size();n++)
    {

      double PINX = solution[n].find(ID_vector[PELTIERNX])->second;
      double PINY = solution[n].find(ID_vector[PELTIERNY])->second;
      double PINZ = solution[n].find(ID_vector[PELTIERNZ])->second;
      double PIPX = solution[n].find(ID_vector[PELTIERPX])->second;
      double PIPY = solution[n].find(ID_vector[PELTIERPY])->second;
      double PIPZ = solution[n].find(ID_vector[PELTIERPZ])->second;
    

      vector_heat_source[n](0) = PINX + PIPX; 
      vector_heat_source[n](1) = PINY + PIPY; 
      vector_heat_source[n](2) = PINZ + PIPZ; 
     
    }

    
  }

 
}






//-------------------------------------------------------------------------//
