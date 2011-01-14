// $Id$


#include "DSSCModel.h"
#include "Database.h"
//#include "Traps.h"


#include "elem.h"


//TIBER_MODULE(DSSCModel, dscbulk, default)


using namespace std;


DSSCModel::DSSCModel(const ModelOptions& options)
  : PhysicalModel(options),
    _porosity(0.5),
    _is_electrolyte(true),
    _is_TiO2(true),
    _electrons("electron"),
    _iodide("iodide"),
    _triiodide("triiodide"),
    _cation("cation"),
    _ke(1e4),
    _beta(1.0),
    _k3(1e8),
    _generation(0.0),
    _alpha(0.2),
    _alpha2(0.0),
    _deltaG(0.0),
    _x0(0.0),
    _perm_ox(85.0),
    _perm_elec(117.0),
    //_permittivity(100.0),
    _elem(NULL)
{
}



void
DSSCModel::do_init(void)
{
  _is_electrolyte = get_option("electrolyte", _is_electrolyte);
  _is_TiO2 = get_option("TiO2", _is_TiO2);

  // prepare porosity for any situation
  get_parameter("porosity", _porosity);
  if (!is_TiO2())
    _porosity = 1.0;
  else if (!is_electrolyte())
    _porosity = 0.0;

  // Maxwell-Garnet model permittivity
  double A1 = _perm_elec + 2*_perm_ox + 2*_porosity*_perm_elec - 2*_porosity*_perm_ox;
  double A2 = _perm_elec + 2*_perm_ox - _porosity*_perm_elec + _porosity*_perm_ox;
  _permittivity = _perm_ox * A1/A2;

  get_parameter("k_e", _ke);
  get_parameter("beta", _beta);
  get_parameter("k_3", _k3);

  get_parameter("permittivity", _permittivity);
  get_parameter("perm_oxide", _perm_ox);
  get_parameter("perm_electrolyte", _perm_elec);

  get_parameter("generation", _generation);
  get_parameter("alpha", _alpha);
  get_parameter("Light", _x0);
  get_parameter("alpha2", _alpha2);
  get_parameter("deltaG", _deltaG);

  get_parameter("ne", _eq_conc.n);
  get_parameter("nI", _eq_conc.I);
  get_parameter("nI3", _eq_conc.I3);

  get_parameter("mu_e", _mobility.n);
  get_parameter("D_I", _mobility.I);
  get_parameter("D_I3", _mobility.I3);
  get_parameter("D_C", _mobility.C);

  _cation.set_particle_charge(1.0);

  //_eq_conc.n = _porosity * get_options().get_option("n_e", _eq_conc.n);
  _eq_conc.n = (1.0 - _porosity) * _eq_conc.n;
  // the following are given in Mol
  double fac = _porosity * Constants::avogadro / 1e3;
  //  _eq_conc.I = get_options().get_option("n_I", _eq_conc.I) * fac;
  //  _eq_conc.I3 = get_options().get_option("n_I3", _eq_conc.I3) * fac;
  _eq_conc.I = _eq_conc.I * fac;
  _eq_conc.I3 = _eq_conc.I3 * fac;
  _eq_conc.C = _eq_conc.I + _eq_conc.I3;

  double kT = Constants::k_B * SimulationOptions::T;
  //  get_parameter("mu_e", _mobility.n);
  //  _mobility.I = get_option("D_I", _mobility.I) / kT;
  //  _mobility.I3 = get_option("D_I3", _mobility.I3) / kT;
  //  _mobility.C = get_option("D_C", _mobility.C) / kT;
  _mobility.I = _mobility.I / kT;
  _mobility.I3 = _mobility.I3 / kT;
  _mobility.C = _mobility.C / kT;

}



void
DSSCModel::read_database(void)
{
  Database& db = get_database();
  db.set_section("recombination_TiO2-electrolyte");
  _ke = db.get("ke",_ke);
  _beta = db.get("beta",_beta);
  _k3 = db.get("k3",_k3);

  db.set_section("permittivity");
  _perm_ox = db.get("permittivity_oxide",_perm_ox);
  _perm_elec = db.get("permittivity_electrolyte",_perm_elec);

  db.set_section("mobility");
  _mobility.n = db.get("mobility_oxide",0.3);
  _mobility.I = db.get("mobility_iodine",8.5e-6);
  _mobility.I3 = db.get("mobility_triiodide",8.5e-6);
  _mobility.C = db.get("mobility_cation",8.5e-6);

  db.set_section("density");
  _eq_conc.n = ("ne",2.3e4);
  _eq_conc.I = ("nI",0.45);
  _eq_conc.I3 = ("nI3",0.05);
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
    double gen1 =  1e4 * _alpha * _generation * exponential;

    double exponential2 = _alpha2 * exp( -_alpha2 * abs(_pd.coordinates(0) - _x0) );
    _pd.generation_rate =  gen1 + _deltaG * exponential2;
    //_pd.generation_rate =  _generation;

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

/*
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
  _pd.recombination_rate_derivatives[1] = -_ke * ( r / _pd.density_I +
     _eq_conc.n * sqrt_I3_I_eq );
  _pd.recombination_rate_derivatives[2] = _ke * r / _pd.density_I3;
*/


/*
  double sqrt_I3_I = sqrt(_pd.density_I3 / _pd.density_I);
  double n_I_p3 = _eq_conc.I *_eq_conc.I * _eq_conc.I;
  double sqrt_I3_I_eq = sqrt(_eq_conc.I3 / n_I_p3);
  // rate
  double r = pow(_pd.density_n,_beta) * sqrt_I3_I;
  double g = pow(_eq_conc.n,_beta) * _pd.density_I * sqrt_I3_I_eq;
  _pd.recombination_rate = _ke * (r - g);

  // derivative
  _pd.recombination_rate_derivatives = vector<double>(4, 0.0);
  _pd.recombination_rate_derivatives[0] = _ke * sqrt_I3_I * _beta * pow(_pd.density_n,_beta-1);
  _pd.recombination_rate_derivatives[1] = -_ke * (0.5 * r / _pd.density_I +
     pow(_eq_conc.n,_beta) * sqrt_I3_I_eq);
  _pd.recombination_rate_derivatives[2] = _ke * 0.5 * r / _pd.density_I3;
*/


  double n0 = _eq_conc.n;
  if (n0 <= _generation/_k3)
  {
     n0 = _generation/_k3;
  }
  if (n0 <= 1e-3)
  {
    n0 = 100;
  }
  double dens_norm = _pd.density_n/n0;
  double dens_norm_dark = _eq_conc.n/n0;

  double dens_beta = pow(_pd.density_n, _beta) * pow(n0, 1.0 - _beta);
  double dens_dark_beta = pow(dens_norm_dark, _beta) * n0;
  double dens_beta_der = _beta * pow(_pd.density_n, _beta - 1) * pow(n0,1-_beta);

  double sqrt_I3_I = sqrt(_pd.density_I3 / _pd.density_I);
  double n_I_p3 = _eq_conc.I *_eq_conc.I * _eq_conc.I;
  double sqrt_I3_I_eq = sqrt(_eq_conc.I3 / n_I_p3);
  // rate
  double r = dens_beta * sqrt_I3_I;
  double g = dens_dark_beta * _pd.density_I * sqrt_I3_I_eq;
  _pd.recombination_rate = _ke * (r - g);
  // derivative
  _pd.recombination_rate_derivatives = vector<double>(4, 0.0);
  _pd.recombination_rate_derivatives[0] = _ke * sqrt_I3_I * dens_beta_der;
  _pd.recombination_rate_derivatives[1] = -_ke * (0.5 * r / _pd.density_I +
     dens_dark_beta * sqrt_I3_I_eq);
  _pd.recombination_rate_derivatives[2] = _ke * 0.5 * r / _pd.density_I3;

}


/*void
DSSCModel::calculate_traps(void)
{
  //double Ec = get_conduction_band_edge() - _pd->electric_potential;
  //double Ev = get_valence_band_edge() - _pd->electric_potential;
  double Ec = 0.93 - _pd->electric_potential;
  double Ev = (0.93 - 3.2)  - _pd->electric_potential;

  _pd->ionized_electron_traps = 0.0;
  _pd->ionized_electron_traps_derivative = 0.0;
  if (_etraps.size() > 0)
  {
    double nt = 0, dnt = 0;
    double kT = Constants::k_B * SimulationOptions::T;
    set<Trap*>::iterator it(_etraps.begin());
    const set<Trap*>::iterator end(_etraps.end());
    for ( ; it != end; ++it)
    {
      (*it)->set_energies(Ec, Ev, -_pd->fermi_n, kT);
      nt += (*it)->get_ionized_density();
      dnt += (*it)->get_ionized_density_derivative();
    }

    _pd->ionized_electron_traps = nt;
    _pd->ionized_electron_traps_derivative = dnt;
  }
}*/

