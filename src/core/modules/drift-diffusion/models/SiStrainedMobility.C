// $Id$

#include "SiStrainedMobility.h"
#include "DriftDiffusionProperties.h"

#include "Database.h"



TIBER_MODULE(SiStrainedMobility, mobility, strained_Si)



SiStrainedMobility::SiStrainedMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    mu0_(1000),
    exp_(1),
    _alpha_c(0), _beta_c(0), _gamma_c(1),
    _alpha_t(0), _beta_t(0), _gamma_t(1)
{
}


void
SiStrainedMobility::read_database(void)
{

  const Database& db = get_database();
  db.set_section("mobility/constant");

  std::vector<double> data(2, 0);
  db.get("mu_max", data, true);
  mu0_ = get_carrier_type() == 'e' ? data[0] : data[1];

  data = std::vector<double>(2, 0);
  db.get("exponent", data);
  exp_ = get_carrier_type() == 'e' ? data[0] : data[1];

}



void
SiStrainedMobility::do_init(void)
{
  get_parameter("mu", mu0_);
  // we allow also mu_e and mu_h
  get_parameter(std::string("mu_") + get_carrier_type(), mu0_);

  // this is for holes only!
  get_parameter("alpha_c", _alpha_c);
  get_parameter("beta_c", _beta_c);
  get_parameter("gamma_c", _gamma_c);
  get_parameter("alpha_t", _alpha_t);
  get_parameter("beta_t", _beta_t);
  get_parameter("gamma_t", _gamma_t);

  std::string sim_name = "";
  get_parameter("strain_simulation", sim_name);
  _strain.set_simulation(sim_name);

}



double
SiStrainedMobility::get_mobility(void)
{
  Tensor2Sym strain(0);
  _strain.get_strain(get_driftdiffusionproperties().get_element(),
      get_driftdiffusionproperties().get_coordinates(), strain);

  // assume always along x direction
  double eps = strain(1,1);

  double mob = mu0_;

  double a = _alpha_c;
  double b = _beta_c;
  double g = _gamma_c;
  if (eps > 0)
  {
    a = _alpha_t;
    b = _beta_t;
    g = _gamma_t;
  }

  double egeps = exp(g*eps);
  mob += a * (1 - egeps) / (1 + b*egeps);


  double T = get_driftdiffusionproperties().get_lattice_temperature();
  return mob * std::pow(T / T0, -exp_);
}



void
SiStrainedMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}


