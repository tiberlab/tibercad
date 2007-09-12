// $Id$

#include "ParticleDensity.h"
#include "SimulationInterface.h"


// GNU scientific library
#include <gsl/gsl_sf_fermi_dirac.h>

using namespace std;




ParticleDensity::ParticleDensity(double particle_charge,
    TiberCad::Statistics statistics)
: _particle_charge(particle_charge),
  _statistics(statistics),
  _use_quantum(false),
  _density_id(INVALID_ID),
  _elem(NULL),
  _density(-1.0),
  _density_derivative(-1.0)
{
}




void
ParticleDensity::add_quantum_density(const std::string& name)
{
  if (name != "")
  {
    SimulationInterface* qd = SimulationInterface::find_simulation(name);
    if (qd == NULL)
    {
      string msg("ParticleDensity: ");
      msg += "no quantum density simulation '" + name + "' found.";
      throw InitFailedException(msg);
    }

    // we assume that the density variable has this name:
    string density_name("density");
    
    _density_id = qd->get_variable_id(density_name);

    if (_density_id == INVALID_ID)
    {
      string msg("ParticleDensity: ");
      msg += "quantum density simulation '" + name +
        "' has no variable '" + density_name + "'";
      throw InitFailedException(msg);
    }

    // at this point we have for sure a quantum density simulation
    _quantum_density.insert(qd);
    use_quantum_density();
  }
}




template <>
inline
void
ParticleDensity::classical_density<TiberCad::BOLTZMANN>(void)
{
  _density = _N_eff * exp(_argument);
}



template <>
inline
void
ParticleDensity::classical_density_derivative<TiberCad::BOLTZMANN>(void)
{
  classical_density<TiberCad::BOLTZMANN>();
  _density_derivative = _density / _kT;
}



template <>
inline
void
ParticleDensity::classical_density<TiberCad::FERMIDIRAC>(void)
{
  _density = _N_eff * gsl_sf_fermi_dirac_half(_argument);
}



template <>
inline
void
ParticleDensity::classical_density_derivative<TiberCad::FERMIDIRAC>(void)
{
  _density_derivative = _N_eff * gsl_sf_fermi_dirac_mhalf(_argument) / _kT;
}



bool
ParticleDensity::quantum_density(void)
{
  set<SimulationInterface*>::iterator it(_quantum_density.begin());
  const set<SimulationInterface*>::iterator end(_quantum_density.end());

  bool flag = false;
  
  _density = 0.0;
  for ( ; it != end; ++it)
    if ((*it)->is_solved())
    {
      flag = true;
      _density += (*it)->get_solution(_elem, _p, _density_id);
    }
    
  return flag;
}


bool
ParticleDensity::quantum_density_derivative(void)
{
  set<SimulationInterface*>::iterator it(_quantum_density.begin());
  const set<SimulationInterface*>::iterator end(_quantum_density.end());

  bool flag = false;
  
  _density_derivative = 0.0;
  for ( ; it != end; ++it)
    if ((*it)->is_solved())
    {
      flag = true;
    }

  return flag;
}



void
ParticleDensity::calculate_density(void)
{
  /* classical density is calculated when
   * * either _use_quantum is false
   * * or quantum_density() returns false
   */
  if (!_use_quantum || !quantum_density())
  {
    switch (_statistics)
    {
      case TiberCad::FERMIDIRAC:
        classical_density<TiberCad::FERMIDIRAC>();
        break;
      default: // Boltzmann
        classical_density<TiberCad::BOLTZMANN>();
        break;
    }
  }
}




void
ParticleDensity::calculate_density_derivative(void)
{
  if (!_use_quantum || !quantum_density_derivative())
  {
    switch (_statistics)
    {
      case TiberCad::FERMIDIRAC:
        classical_density_derivative<TiberCad::FERMIDIRAC>();
        break;
      default: // Boltzmann
        classical_density_derivative<TiberCad::BOLTZMANN>();
        break;
    }
  }
}


