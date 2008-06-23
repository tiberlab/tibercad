// $Id$

#include "TestMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"


TIBER_MODULE(TestMobility, test)


void
TestMobility::read_database(void)
{
}



void
TestMobility::do_init(void)
{
  _mu_min = get_parameter("mu_min", _mu_min);
  _mu_max = get_parameter("mu_max", _mu_max);
  _E_max = get_parameter("E_max", _E_max);

  std::string force = get_parameter("driving_force", "grad_fermi");
  if (force == "efield")
    _force = EFIELD;
  else if (force == "grad_fermi")
    _force = GRADFERMI;
  else
  {
    std::string msg("TestMobility: Unknown driving force '");
    msg += force + "'.";
    throw InitFailedException(msg);
  }
    
}



double
TestMobility::get_mobility(void)
{
  double E = 0.0;
  if (get_carrier_type() == 'e')
    E = get_driftdiffusionproperties().get_grad_fermi_e().size();
  else
    E = get_driftdiffusionproperties().get_grad_fermi_h().size();
 
  double mu = _mu_max + (_mu_min - _mu_max) / _E_max * std::min(E, _E_max);
  return mu;
}



void
TestMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}


void
TestMobility::get_derivative_grad_fermi(RealGradient& dm)
{
  dm.zero();

  double E = 0.0;
  if (get_carrier_type() == 'e')
    E = get_driftdiffusionproperties().get_grad_fermi_e().size();
  else
    E = get_driftdiffusionproperties().get_grad_fermi_h().size();
  
  if ((_force == GRADFERMI) && (E > 1.0))
  {

    double dmu = 0;

    if (E < _E_max)
    {
      dmu = (_mu_min - _mu_max) / _E_max / E;
      
      if (get_carrier_type() == 'e')
        dm = dmu * get_driftdiffusionproperties().get_grad_fermi_e();
      else
        dm = dmu * get_driftdiffusionproperties().get_grad_fermi_h();
    }

  }
}


void
TestMobility::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
}

