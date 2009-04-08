#include "EigenvalueProblem.h"
#include "Constants.h"

void EigenvalueProblem::get_eigenvalues(std::vector<double>& values) const
{

  unsigned int n = _solution.size();
  values.resize(n);

  for (unsigned int i = 0; i < n; i++)
  {
    values[i] = _solution[i].eigen_energy;  
  }  

}

void EigenvalueProblem::get_populations(std::vector<double>& values) const
{
 
  unsigned int n = _solution.size();
  values.resize(n);
 
  for (unsigned int i = 0; i < n; i++)
  {
    if(_solution[i].statistics == "Fermi")
    {
      values[i] = Fermi(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].Temperature);

      if(_solution[i].particle == "hl" || _solution[i].particle == "hole")
      {
	values[i] = 1 - values[i];
      }
  
    }
    else
    {
      values[i] = Bose(_solution[i].eigen_energy, _solution[i].electro_chem_pot, 
		       _solution[i].Temperature);
      
    }
  }
 
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
