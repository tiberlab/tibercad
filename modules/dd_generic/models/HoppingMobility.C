// $Id: HoppingMobility.C $

#include "HoppingMobility.h"
#include "DriftDiffusionProperties.h"

#include "tibercad/module/TiberModule.h"

#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"


void
HoppingMobility::do_init(void)
{
  _sigma = get_option("sigma", _sigma);
  _N0 = get_option("N0", _N0);
  _nu0 = get_option("nu0", _nu0);
}

double HoppingMobility::get_mobility(void)
{
  ID id = get_carrier();
  double dens = get_driftdiffusionproperties().get_q_density(id);
  
  double F = get_driftdiffusionproperties().get_electric_field().norm();
  double kT = get_driftdiffusionproperties().get_lattice_temperature();

  const double c1 = 1.8e-9;
  const double c2 = 0.42;
  	
  double s = _sigma / kT;
  double a = pow(_N0, (-1.0/3.0));
  
  double mu0 = _nu0 * a * a / _sigma;
    
  double mu0s = mu0 * c1 * exp(-1.0 * c2 * s * s);

  double n_dens = dens / _N0;

  
  if (n_dens > 0.1)  //cfr Phys. Rev. B 78 (2008) 085207 for constraints on n_dens and n_F  //0.1
  {
    n_dens = 0.1;  //0.1
  }
  

  double delta = 2 * (log(s*s - s) - log(log(4))) / (s*s);
  double mu0s_n = mu0s * exp(0.5 * (s*s - s) * pow(2* n_dens, delta));

  double n_F = F * a / _sigma;

  
  if (n_F > 2)  // 2
  {
    n_F = 2;
  }
  

  double fsF = exp(0.44 * (sqrt(s*s*s) - 2.2) * (sqrt(1 + 0.8*n_F*n_F) - 1));
  //std::cout<<"s="<<s<<" n_dens="<<n_dens<<" n_F="<<n_F<<" mu="<<mu0s_n / mu0<<std::endl;
  return mu0s_n * fsF;
}

double
HoppingMobility::get_derivative_potential(void)
{
  double kT = get_driftdiffusionproperties().get_lattice_temperature();
  double mu, dens, dens_der;

  ID id = get_carrier();
  mu = get_driftdiffusionproperties().get_q_mobility(id);
  dens = get_driftdiffusionproperties().get_q_density(id);
  dens_der = get_driftdiffusionproperties().get_q_density_derivative(id);

  double n_dens = dens / _N0;
  if (n_dens > 0.1)
  {
    return 0.0;
  }
  else
  {
    double s = _sigma / kT;
    double delta = 2 * (log(s*s - s) - log(log(4))) / (s*s);

    return mu * 0.5 * delta * (s*s - s) * pow(2 * n_dens, delta) * dens_der / dens;
  }
}

void
HoppingMobility::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  ID id = get_carrier();
  double mu = this->get_mobility();

  double F = get_driftdiffusionproperties().get_electric_field().norm();
  double kT = get_driftdiffusionproperties().get_lattice_temperature();
  double s = _sigma / kT;
  double a = pow(_N0, (-1.0/3.0));

  double n_F = F * a / _sigma;
  double dmu;
  if (n_F > 2)
  {
    dm.zero();
  }
  else
  {
    dmu = 2 * mu * 0.22 * (sqrt(s*s*s) -2.2) * 0.8 * a * a / (_sigma * _sigma * sqrt(1 + 0.8 * n_F * n_F));
    dm = -1.0 * dmu * get_driftdiffusionproperties().get_electric_field();
  }
}


void
HoppingMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
}
