// $Id$


#include "DSSCModel.h"


#include "elem.h"


using namespace std;


DSSCModel::DSSCModel(void)
  : _porosity(0.5),
    _is_electrolyte(true),
    _is_TiO2(true),
    _electrons("electron"),
    _iodide("iodide"),
    _triiodide("triiodide"),
    _cation("cation"),
    _ke(0.0),
    _k3(1.0),
    _generation(0.0),
    _alpha(0.0),
    _x0(0.0),
    _alpha2(0.0),
    _deltaG(0.0)
{
}



void
DSSCModel::do_init(void)
{
  _is_electrolyte = get_option("electrolyte", _is_electrolyte);
  _is_TiO2 = get_option("TiO2", _is_TiO2);

  // prepare porosity for any situation
  get_parameter("porosity", _porosity);
  //if (!is_TiO2())
  //  _porosity = 0.0;
  //else if (!is_electrolyte())
  //  _porosity = 1.0;

  _cation.set_particle_charge(1.0);

  _eq_conc.n = _porosity * get_options().get_option("n_e", _eq_conc.n);
  // the following are given in Mol
  double fac = (1.0 - _porosity) * Constants::avogadro / 1e3;
  _eq_conc.I = get_options().get_option("n_I", _eq_conc.I) * fac;
  _eq_conc.I3 = get_options().get_option("n_I3", _eq_conc.I3) * fac;
  _eq_conc.C = _eq_conc.I + _eq_conc.I3;

  double kT = Constants::k_B * SimulationOptions::T;
  get_parameter("mu_e", _mobility.n);
  _mobility.I = get_option("D_I", _mobility.I) / kT;
  _mobility.I3 = get_option("D_I3", _mobility.I3) / kT;
  _mobility.C = get_option("D_C", _mobility.C) / kT;

  get_parameter("k_e", _ke);
  get_parameter("k_3", _k3);

  get_parameter("permittivity", _permittivity);

  get_parameter("generation", _generation);

  get_parameter("alpha", _alpha);

  get_parameter("alpha2", _alpha2);
  get_parameter("deltaG", _deltaG);
}



void
DSSCModel::reinit(const Elem* elem)
{
  if (_elem != elem)
  {
    _elem = elem;
    _pd.coordinates = elem->centroid();
    _pd.kT = Constants::k_B * _lattice_temp.get_temperature(elem, _pd.coordinates);

    prepare_element_data();
  }
}


void
DSSCModel::copy_from(const PhysicalModelInterface* rhs)
{
}


void
DSSCModel::do_print_info(void)
{
  cerr << "n_e_0 = " << _eq_conc.n << "  mu_n = " << _mobility.n << endl;
  cerr << "n_I_0 = " << _eq_conc.I << "  mu_I = " << _mobility.I << endl;
  cerr << "n_I3_0 = " << _eq_conc.I3 << "  mu_I3 = " << _mobility.I3 << endl;
  cerr << "n_C_0 = " << _eq_conc.C << "  mu_C = " << _mobility.C << endl;
}



void
DSSCModel::calculate_densities(void)
{
  _pd.density_n = 0.0;
  _pd.generation_rate = 0.0;
  if (is_TiO2())
  {
    _electrons.set_element_and_point(_elem, _pd.coordinates);
    _electrons.set_classical_parameters(_eq_conc.n,
        -_pd.electric_potential, -_pd.fermi_n, _pd.kT);
    _pd.density_n = _electrons.get_particle_density();

    // generation has to be calculated here
    double exponential = exp( -_alpha * abs(_pd.coordinates(0) - _x0) );
    double gen1 =  _generation * exponential;
    //if (_generation > _deltaG)
    //{
       double exponential2 = _alpha2 * exp( -_alpha2 * abs(_pd.coordinates(0) - _x0) );
       _pd.generation_rate =  gen1 + _deltaG * exponential2;
    //}
    //else
    //{
    //   _pd.generation_rate =  gen1;
    //}

  }

  _pd.density_I = _pd.density_I3 = _pd.density_C = 0.0;
  if (is_electrolyte())
  {
    _iodide.set_element_and_point(_elem, _pd.coordinates);
    _iodide.set_classical_parameters(_eq_conc.I, -_pd.electric_potential,
        -_pd.fermi_I, _pd.kT);
    _pd.density_I = _iodide.get_particle_density();

    _triiodide.set_element_and_point(_elem, _pd.coordinates);
    _triiodide.set_classical_parameters(_eq_conc.I3, -_pd.electric_potential,
       -_pd.fermi_I3, _pd.kT);
    _pd.density_I3 = _triiodide.get_particle_density();

    _cation.set_element_and_point(_elem, _pd.coordinates);
    _cation.set_classical_parameters(_eq_conc.C, _pd.electric_potential,
        _pd.fermi_C, _pd.kT);
    _pd.density_C = _cation.get_particle_density();
  }

  _pd.ionized_dye = _pd.generation_rate / _k3;
}



void
DSSCModel::calculate_net_recombination_rate(void)
{
  double sqrt_I3_I = sqrt(_pd.density_I3 / _pd.density_I);
  double n_I_p3 = _eq_conc.I *_eq_conc.I * _eq_conc.I;
  double sqrt_I3_I_eq = sqrt(_eq_conc.I3 / n_I_p3);
  // rate
  double r = _pd.density_n * sqrt_I3_I;
  double g = _eq_conc.n * _pd.density_I * sqrt_I3_I_eq;
  _pd.recombination_rate = _ke * (r - g);

  // derivative
  _pd.recombination_rate_derivatives = vector<double>(4, 0.0);
  _pd.recombination_rate_derivatives[0] = _ke * sqrt_I3_I;
  _pd.recombination_rate_derivatives[1] = -_ke * (0.5 * r / _pd.density_I +
     _eq_conc.n * sqrt_I3_I_eq);
  _pd.recombination_rate_derivatives[2] = _ke * 0.5 * r / _pd.density_I3;
}
