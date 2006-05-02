#include "StrainedSemiconductorModel.h"
#include "DDsemiconductor.h"

#include "elem.h"
#include "macrostrain.h"
#include "tensor.h"

#include <iostream>

using namespace DriftDiffusionDefs;

StrainedSemiconductorModel::StrainedSemiconductorModel(
    Macrostrain* strain)
  : SemiconductorModel(),
    _strain(strain),
    _ignore_strain(false)
{
}


StrainedSemiconductorModel::StrainedSemiconductorModel(
    const StrainedSemiconductorModel& model)
  : SemiconductorModel(model),
    _strain(model._strain),
    _ignore_strain(model._ignore_strain)
{
}

void
StrainedSemiconductorModel::prepare_element_data(void)
{
  if (!_ignore_strain)
  {
    assert(elem != NULL);
    if ((elem != NULL))
    {
      Tensor1 pol = _strain->get_piezopolarization(elem);
      polarization(0) = pol(1);
      polarization(1) = pol(2);
      polarization(2) = pol(3);

      // set strain
      //get_physical_model()->set_strain(Tensor2Sym(0));
      get_physical_model()->set_strain(_strain->get_strain(elem, true));
      cout << "strain: " << norm(_strain->get_strain(elem, true)) << endl;
      cout << setw(15) << _strain->get_strain(elem, true) << endl;
      cout << setw(15) << _strain->get_strain(elem) << endl;

      // calculate new equilibrium values
      get_physical_model()->calculate_conduction_band_extremum();
      get_physical_model()->calculate_valence_band_extremum();

      // get all band properties
      extract_band_properties();

      // this sets the band edges and the effective DOS in the base class
      SemiconductorModel::prepare_element_data();

      // calculate the new equlibrium electro-chemical potential
      double Ec = get_conduction_band_edge();
      double Eg = Ec - get_valence_band_edge();
      double Ef = Ec - Ec0 + Ef0 - (Eg - Eg0) / gamma;
      equilibrium_fermi_level = Ef;

      // calculate equilibrium electron and hole concentration
      // the electrons:
      double kT = electron_vt;
      double arg = (Ef - Ec) / kT;
      double n = get_conduction_band_properties().effective_DOS
        * density<TiberCad::FERMIDIRAC>(arg);
      equilibrium_electron_density = n;
      // the holes:
      arg = (Ec - Eg - Ef) / kT;
      double p = get_valence_band_properties().effective_DOS
        * density<TiberCad::FERMIDIRAC>(arg);
      equilibrium_hole_density = p;
    }
  }
  else
    SemiconductorModel::prepare_element_data();
}

void
StrainedSemiconductorModel::calculate_equilibrium_properties(
    int coupling, double temperature)
{

  // explicitly clear strain
  //get_physical_model()->set_strain(Tensor2Sym(0));

  // call method of parent class
  SemiconductorModel::calculate_equilibrium_properties(coupling, temperature);

  Ef0 = get_equilibrium_fermi_level();
  Ec0 = get_conduction_band_edge();
  Eg0 = Ec0 - get_valence_band_edge();

  /*
   *              d/dEf (n - Nd+)
   * gamma = 1 + -----------------
   *              d/dEf (-p + Na+)
   */
  double a = get_electron_density_derivative()
    - get_ionized_donor_density_derivative();
  double b = get_ionized_acceptor_density_derivative()
    - get_hole_density_derivative();
  
  // b can go to zero, e.g. if Na << Nd
  if (b < 1e-12 * a)
    gamma = 1e12;
  else
    gamma = 1 + a / b;

}
