// $Id$

#include "ConstantMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"


#ifdef TIBER_MODULE

# undef TIBER_MODELNAME
# define TIBER_MODELNAME constant

extern "C"
{

  ConstantMobility* create(void)
  {
    return new ConstantMobility();
  }

  void destroy(PhysicalModelInterface* p)
  {
    delete p;
  }
}

#endif




void
ConstantMobility::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());
  
  std::string s("mu_max_");
  s += get_carrier_type();
  _mu0 = data(s.c_str(), _mu0);

  s = "exponent_";
  s += get_carrier_type();
  _exp = data(s.c_str(), _exp);
}



void
ConstantMobility::do_init(void)
{
  _mu0 = get_options().get_option("mu", _mu0);
  std::string s("mu_");
  s += get_carrier_type();
  _mu0 = get_material()->get_options().get_option(s, _mu0);
}



double
ConstantMobility::get_mobility(void)
{
  double T = get_driftdiffusionproperties().get_lattice_temperature();
  return _mu0 * std::pow(T / T0, -_exp);
}



void
ConstantMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}



void
ConstantMobility::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const ConstantMobility* scA =
    dynamic_cast<const ConstantMobility*>(comp_A);
  const ConstantMobility* scB =
    dynamic_cast<const ConstantMobility*>(comp_B);

  _mu0 = alloy(scA->_mu0, scB->_mu0, xa);
  _exp = alloy(scA->_exp, scB->_exp, xa);
}

