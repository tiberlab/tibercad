/*  
 * This file is part of the tiberCAD module driftdiffusion.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file HoppingMobility.C
 * \brief tiberCAD driftdiffusion module implementation.
 *
 * \note This file is part of module driftdiffusion.
 */


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
  double dens;
  if (get_carrier_type() == 'e')
  {
    dens = get_driftdiffusionproperties().get_electron_density();
  }
  else
  {
    dens = get_driftdiffusionproperties().get_hole_density();
  }
  
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
  if (get_carrier_type() == 'e')
  {
    mu = get_driftdiffusionproperties().get_electron_mobility();
    dens = get_driftdiffusionproperties().get_electron_density();
    dens_der = get_driftdiffusionproperties().get_electron_density_derivative();
  }
  else
  {
    mu = get_driftdiffusionproperties().get_hole_mobility();
    dens = get_driftdiffusionproperties().get_hole_density();
    dens_der = get_driftdiffusionproperties().get_hole_density_derivative();
  }
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
  double mu;
  if (get_carrier_type() == 'e')
  {
    mu = get_driftdiffusionproperties().get_electron_mobility();
  }
  else
  {
    mu = get_driftdiffusionproperties().get_hole_mobility();
  }
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
HoppingMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}

void
HoppingMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
}
