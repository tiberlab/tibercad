// $Id: HGCoupling.C 3542 2013-03-01 09:31:59Z maufder $

#include "HGCoupling.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"


using namespace std;




void
HGCoupling::do_init(void)
{
  get_parameter("C", C_);

}



void
HGCoupling::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  //*
  if (coupling() == "eHG")
  {
    double n  = dd.get_electron_density();
    double nG  = dd.get_electron_density_guest();
    double Efn = dd.get_electron_electro_chemical_potential();
    double EfnG = dd.get_electron_electro_chemical_potential_guest();
    double Ec = dd.get_conduction_band_edge();
    double EcG = dd.get_conduction_band_guest().get_band_edge();
    double kT = dd.get_lattice_temperature();

    double arg = (Ec - EcG) / kT;
    //double thermal_activation1 = (arg > 0) ? exp(-arg) : 1;
    double thermal_activation = (arg < 0) ? exp(arg) : 1;
    //thermal_activation2 = 1;

    double fac = (1 - exp((Efn - EfnG) / kT));
    //double fac = fabs(tanh((Efn - EfnG) / kT));
    recomb_e = C_ * (n * (1e18 - nG) ) * fac * thermal_activation;
  }
  else
  {
    double p  = dd.get_hole_density();
    double pG  = dd.get_hole_density_guest();
    double Efp = dd.get_hole_electro_chemical_potential();
    double EfpG = dd.get_hole_electro_chemical_potential_guest();
    double Ev = dd.get_valence_band_edge();
    double EvG = dd.get_valence_band_guest().get_band_edge();
    double kT = dd.get_lattice_temperature();

    double arg = (EvG - Ev) / kT;
    //double thermal_activation1 = (arg < 0) ? exp(arg) : 1;
    double thermal_activation = (arg < 0) ? exp(arg) : 1;
    //thermal_activation2 = 1;

    double fac = (1 - exp((Efp - EfpG) / kT));
    //double fac = fabs(tanh((Efp - EfpG) / kT));
    recomb_e = -C_ * (p * (1e18 - pG) ) * fac * thermal_activation;
  }
  //*/
  /*
  if (coupling() == "eHG")
  {
    double n  = dd.get_electron_density();
    double nG  = dd.get_electron_density_guest();
    double Efn = dd.get_electron_electro_chemical_potential();
    double EfnG = dd.get_electron_electro_chemical_potential_guest();
    double Ec = dd.get_conduction_band_edge();
    double EcG = dd.get_conduction_band_guest().get_band_edge();
    double kT = dd.get_lattice_temperature();

    double arg = (EcG - Ec) / kT;
    double thermal_activation1 = (arg > 0) ? exp(-arg) : 1;
    double thermal_activation2 = (arg < 0) ? exp(arg) : 1;
    thermal_activation2 = 1;

    //double fac = (1 - exp((Efn - EfnG) / kT));
    double fac = fabs(tanh((Efn - EfnG) / kT));
    recomb_e = C_ * (n * (1e18 - nG) * thermal_activation1 - nG * 1e19 * thermal_activation2) * fac;
  }
  else
  {
    double p  = dd.get_hole_density();
    double pG  = dd.get_hole_density_guest();
    double Efp = dd.get_hole_electro_chemical_potential();
    double EfpG = dd.get_hole_electro_chemical_potential_guest();
    double Ev = dd.get_valence_band_edge();
    double EvG = dd.get_valence_band_guest().get_band_edge();
    double kT = dd.get_lattice_temperature();

    double arg = (EvG - Ev) / kT;
    double thermal_activation1 = (arg < 0) ? exp(arg) : 1;
    double thermal_activation2 = (arg > 0) ? exp(-arg) : 1;
    thermal_activation2 = 1;

    //double fac = (1 - exp((Efn - EfnG) / kT));
    double fac = fabs(tanh((Efp - EfpG) / kT));
    recomb_e = -C_ * (p * (1e18 - pG) * thermal_activation1 - pG * 1e19 * thermal_activation2) * fac;
  }
  */
  recomb_h = - recomb_e;
}



void
HGCoupling::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double kT = dd.get_lattice_temperature();
  //*
  if (coupling() == "eHG")
  {
    double n  = dd.get_electron_density();
    double nG  = dd.get_electron_density_guest();
    double Efn = dd.get_electron_electro_chemical_potential();
    double EfnG = dd.get_electron_electro_chemical_potential_guest();
    double Ec = dd.get_conduction_band_edge();
    double EcG = dd.get_conduction_band_guest().get_band_edge();

    double arg = (Ec - EcG) / kT;
    //double thermal_activation1 = (arg > 0) ? exp(-arg) : 1;
    double thermal_activation = (arg < 0) ? exp(arg) : 1;
    //thermal_activation2 = 1;

    double fac = (1 - exp((Efn - EfnG) / kT));
    //double fac = fabs(tanh((Efn - EfnG) / kT));

    recomb_e[0] =  C_ * (1e18 - nG) * fac * thermal_activation; // dR/dn
    recomb_e[1] = -C_ * n * fac * thermal_activation; // dR/dnG
    recomb_e[2] = -exp((Efn - EfnG) / kT) / kT * C_ * (n * (1e18 - nG) ) * thermal_activation;
    recomb_e[3] =  exp((Efn - EfnG) / kT) / kT * C_ * (n * (1e18 - nG) ) * thermal_activation;
    recomb_h[0] = -C_ * (1e18 - nG) * fac * thermal_activation; // dR/dn
    recomb_h[1] =  C_ * n * fac * thermal_activation; // dR/dnG
    recomb_h[2] =  exp((Efn - EfnG) / kT) / kT * C_ * (n * (1e18 - nG) ) * thermal_activation;
    recomb_h[3] = -exp((Efn - EfnG) / kT) / kT * C_ * (n * (1e18 - nG) ) * thermal_activation;
  }
  else
  {
    double p  = dd.get_hole_density();
    double pG  = dd.get_hole_density_guest();
    double Efp = dd.get_hole_electro_chemical_potential();
    double EfpG = dd.get_hole_electro_chemical_potential_guest();
    double Ev = dd.get_valence_band_edge();
    double EvG = dd.get_valence_band_guest().get_band_edge();

    double arg = (EvG - Ev) / kT;
    //double thermal_activation1 = (arg < 0) ? exp(arg) : 1;
    double thermal_activation = (arg < 0) ? exp(arg) : 1;
    //thermal_activation2 = 1;

    //double fac = (1 - exp((Efn - EfnG) / kT));
    double fac = (1 - exp((Efp - EfpG) / kT));

    recomb_e[0] = -C_ * (1e18 - pG) * fac * thermal_activation; // dR/dp
    recomb_e[1] =  C_ * p * fac * thermal_activation; // dR/dpG
    recomb_e[2] =  exp((Efp - EfpG) / kT) / kT * C_ * (p * (1e18 - pG) ) * thermal_activation;
    recomb_e[3] = -exp((Efp - EfpG) / kT) / kT * C_ * (p * (1e18 - pG) ) * thermal_activation;
    recomb_h[0] =  C_ * (1e18 - pG) * fac * thermal_activation; // dR/dp
    recomb_h[1] = -C_ * p * fac * thermal_activation; // dR/dpG
    recomb_h[2] = -exp((Efp - EfpG) / kT) / kT * C_ * (p * (1e18 - pG) ) * thermal_activation;
    recomb_h[3] =  exp((Efp - EfpG) / kT) / kT * C_ * (p * (1e18 - pG) ) * thermal_activation;
  }
  //*/
  /*
  if (coupling() == "eHG")
  {
    double n  = dd.get_electron_density();
    double nG  = dd.get_electron_density_guest();
    double Efn = dd.get_electron_electro_chemical_potential();
    double EfnG = dd.get_electron_electro_chemical_potential_guest();
    double Ec = dd.get_conduction_band_edge();
    double EcG = dd.get_conduction_band_guest().get_band_edge();

    double arg = (EcG - Ec) / kT;
    double thermal_activation1 = (arg > 0) ? exp(-arg) : 1;
    double thermal_activation2 = (arg < 0) ? exp(arg) : 1;
    thermal_activation2 = 1;

    //double fac = (1 - exp((Efn - EfnG) / kT));
    double fac = fabs(tanh((Efn - EfnG) / kT));

    recomb_e[0] = C_ * (1e18 - nG) * thermal_activation1 * fac; // dR/dn
    recomb_e[1] = -C_ * (n * thermal_activation1 + 1e19* thermal_activation2) * fac; // dR/dnG
    recomb_e[2] = 0; -C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
    recomb_e[3] = 0; C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
    recomb_h[0] = -C_ * (1e19 - nG) * thermal_activation1 * fac; // dR/dn
    recomb_h[1] = C_ * (n * thermal_activation1 + 1e19 * thermal_activation2) * fac; // dR/dnG
    //recomb_h[0] = -C_ * (1e19 - nG) * fac; // dR/dn
    //recomb_h[1] = C_ * n * fac; // dR/dnG
    recomb_h[2] = 0; C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
    recomb_h[3] = 0; -C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
  }
  else
  {
    double n  = dd.get_hole_density();
    double nG  = dd.get_hole_density_guest();
    double Efn = dd.get_hole_electro_chemical_potential();
    double EfnG = dd.get_hole_electro_chemical_potential_guest();
    double Ev = dd.get_valence_band_edge();
    double EvG = dd.get_valence_band_guest().get_band_edge();

    double arg = (EvG - Ev) / kT;
    double thermal_activation1 = (arg < 0) ? exp(arg) : 1;
    double thermal_activation2 = (arg > 0) ? exp(-arg) : 1;
    thermal_activation2 = 1;

    //double fac = (1 - exp((Efn - EfnG) / kT));
    double fac = -fabs(tanh((Efn - EfnG) / kT));

    recomb_e[0] = C_ * (1e18 - nG) * thermal_activation1 * fac; // dR/dn
    recomb_e[1] = -C_ * (n * thermal_activation1 + 1e19 * thermal_activation2) * fac; // dR/dnG
    recomb_e[2] = 0; -C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
    recomb_e[3] = 0; C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
    recomb_h[0] = -C_ * (1e18 - nG) * thermal_activation1 * fac; // dR/dn
    recomb_h[1] = C_ * (n * thermal_activation1 + 1e19 * thermal_activation2) * fac; // dR/dnG
    //recomb_h[0] = -C_ * (1e19 - nG) * fac; // dR/dn
    //recomb_h[1] = C_ * n * fac; // dR/dnG
    recomb_h[2] = 0; C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
    recomb_h[3] = 0; -C_ / kT * n * (1e18 - nG) * exp((Efn - EfnG) / kT);
  }
  */
}


