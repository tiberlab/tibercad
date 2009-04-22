#include "EigenvalueProblem.h"
#include "Constants.h"

void EigenvalueProblem::get_eigenvalues(const std::string& particle, 
					std::vector<double>& values) const
{

  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);

  for (unsigned int i = 0; i < n; i++)
  {
    if(_solution[i].particle == particle)
    {  
      num_st++;

      values.push_back( _solution[i].eigen_energy ); 
    }
  }  

  values.resize(num_st);

}

unsigned int EigenvalueProblem::get_num_states(const std::string& particle)
{
  unsigned int num_i_states = 0;
  for(unsigned int i=0; i<_solution.size(); i++)
  {
    if(_solution[i].particle == particle) num_i_states++;  
  }
  
  return num_i_states;
}

void EigenvalueProblem::get_populations(const std::string& particle, 
					std::vector<double>& values) const
{
 
  unsigned int n = _solution.size();
  unsigned int num_st = 0;
  values.reserve(n);
 
  for (unsigned int i = 0; i < n; i++)
  {
    if(_solution[i].particle == particle)
    {  
      num_st++;

      if(_solution[i].statistics == "Fermi")
      {      
	double val = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].temperature);

	if(particle == "el" || particle == "electron")
	{
	  values.push_back(val);	  
	}	
	
	if(particle == "hl" || particle == "hole")
	{
	  values.push_back(1-val);	  
	}

      }
      else
      {
	double val = Bose(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].temperature);

	values.push_back(val);	
	
      }

    }
     
  }

  values.resize(num_st);
 
} 


double  EigenvalueProblem::Fermi(double Energy, double Fermi_energy, double Temperature) const
{
  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double el_fermi;
  
  if (exp_arg > 35) 
    el_fermi = 0.0;
  else
    el_fermi = 1.0/(std::exp(exp_arg) + 1.0);
  
  return el_fermi;

}

double  EigenvalueProblem::Bose(double Energy, double electro_chem_pot, double Temperature) const
{
  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - electro_chem_pot)/T_EV;
  
  double bose;
  
  if (exp_arg > 35) 
    bose = 0.0;
  else
    bose = 1.0/(std::exp(exp_arg) - 1.0);
  
  return bose;

}


void EigenvalueProblem::write_states(void)
{

  int num_st=_solution.size();

  for(int i=0; i< num_st; i++)
  {
    std::cerr<<_solution[i].particle<<" "<<_solution[i].eigen_energy<<" "
	     <<_solution[i].statistics<<" "<<_solution[i].electro_chem_pot
	     <<std::endl;
  }


}
