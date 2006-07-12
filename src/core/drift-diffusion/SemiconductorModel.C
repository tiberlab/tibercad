// $Id$

#include "Alloy.h"

#include "SemiconductorModel.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "getpot.h"

#include <iostream>

using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  delete _bulk_model;
}

SemiconductorModel::SemiconductorModel(void)
  : Parent(),
    _bulk_model(NULL),
    _is_prepared(false),
    _exciton_gen_param(0.0)
{
}

/*
SemiconductorModel::SemiconductorModel(
    const SemiconductorModel& model)
  : Parent(model),
    _filename(model._filename),
    _is_prepared(false),
    _exciton_gen_param(model._exciton_gen_param)
{
}
*/

void
SemiconductorModel::read_database(const Dummy&)
{
  GetPot data(_filename);

  const std::string structure = data("structure", "zb");

  if (structure == "zb")
  {
    // defaults for GaAs
    ZbDDsemiconductor::ZbDDparameters params;

    params.EgGamma = data("Eg_G", 1.519);
    params.EgL = data("Eg_L", 1.815);
    params.EgX = data("Eg_X", 1.981);
    params.Ev = data("E_v", 1.346);

    params.m_G = data("m_G", 0.067);
    params.m_l_L = data("m_L_l", 1.9);
    params.m_t_L = data("m_L_t", 0.0754);
    params.m_l_X = data("m_X_l", 1.3);
    params.m_t_X = data("m_X_t", 0.23);

    params.a_c = data("a_c", -9.36);
    params.a_v = data("a_v", -1.21);
    params.b = data("b", -2.0);
    params.d = data("d", -4.8);

    params.delta = data("delta", 0.341);
    params.gamma1 = data("gamma1", 6.98);
    params.gamma2 = data("gamma2", 2.06);
    params.gamma3 = data("gamma3", 2.93);

    params.def_vol_X = data("abs_def_pot_X", -0.16);
    params.def_uniax_X = data("uniax_def_pot_X", 14.26);
    params.def_vol_L = data("abs_def_pot_L", -4.91);
    params.def_uniax_L = data("uniax_def_pot_L", 6.5);

    ZbDDsemiconductor* zbsc = new ZbDDsemiconductor(params);
    _bulk_model = zbsc;
    zbsc->energy_cutoff = 4.0;

    permittivity = data("permittivity", 12.93);
    _e_mobility = data("electron_mobility", 1000.0);
    _h_mobility = data("hole_mobility", 200.0);

  }
  else
  {
    // defaults for GaN
    WzDDsemiconductor::WzDDparameters params;

    params.EgGamma = data("Eg_G", 3.51);
    params.Ev = data("E_v", -0.726);

    params.m_c_zz = data("m_c_zz", 0.20);
    params.m_c_xx = data("m_c_xx", 0.20);

    params.A1 = data("A1", -7.21);
    params.A2 = data("A2", -0.44);
    params.A3 = data("A3", 6.68);
    params.A4 = data("A4", -3.46);
    params.A5 = data("A5", -3.40);
    params.A6 = data("A6", -4.90); 

    params.a_x = data("a_x", -4.9);
    params.a_z = data("a_z", -11.3);

    params.D1 = data("D1", -3.7);
    params.D2 = data("D2", 4.5);
    params.D3 = data("D3", 8.2);
    params.D4 = data("D4", -4.1);
    params.D5 = data("D5", -4.0);
    params.D6 = data("D6", -5.5);
    params.delta_s = data("delta_s", 0.017);
    params.delta_cr = data("delta_cr", 0.010);

    WzDDsemiconductor* wzsc = new WzDDsemiconductor(params);
    _bulk_model = wzsc;
    wzsc->energy_cutoff = 4.0;

    permittivity = data("permittivity", 9.5);
    _e_mobility = data("electron_mobility", 1000.0);
    _h_mobility = data("hole_mobility", 200.0);
  }
}

void
SemiconductorModel::build_alloy(const std::string& component2,
    const std::string& bowing_params, double content)
{
  GetPot data(component2);
  GetPot bowing(bowing_params);

  const std::string structure = data("structure", "zb");

  double (*alloy)(double, double, double, double) =
    Alloy::calculate_VCA_parameter;

  if (structure == "zb")
  {
    ZbDDsemiconductor* sc = dynamic_cast<ZbDDsemiconductor*>(_bulk_model);
    ZbDDsemiconductor::ZbDDparameters& params = sc->get_parameters();

    params.EgGamma = alloy(data("Eg_G", 1.519), params.EgGamma, content,
       bowing("Eg_G", 0.0));
    params.EgL = alloy(data("Eg_L", 1.815), params.EgL, content,
       bowing("Eg_L", 0.0));
    params.EgX = alloy(data("Eg_X", 1.981), params.EgX, content,
       bowing("Eg_X", 0.0));
    params.Ev = alloy(data("E_v", 1.346), params.Ev, content,
       bowing("E_v", 0.0));

    params.m_G = alloy(data("m_G", 0.067), params.m_G, content,
       bowing("m_G", 0.0));
    params.m_t_L = alloy(data("m_L_t", 0.0754), params.m_t_L, content,
       bowing("m_L_t", 0.0));
    params.m_l_L = alloy(data("m_L_l", 1.9), params.m_l_L, content,
       bowing("m_L_l", 0.0));
    params.m_t_X = alloy(data("m_X_t", 1.3), params.m_t_X, content,
       bowing("m_X_t", 0.0));
    params.m_l_X = alloy(data("m_X_l", 0.23), params.m_l_X, content,
       bowing("m_X_l", 0.0));

    params.a_c = alloy(data("a_c", -9.36), params.a_c, content,
       bowing("a_c", 0.0));
    params.a_v = alloy(data("a_v", -1.21), params.a_v, content,
       bowing("a_v", 0.0));
    params.b = alloy(data("b", -2.0), params.b, content,
       bowing("b", 0.0));
    params.d = alloy(data("d", -4.8), params.d, content,
       bowing("d", 0.0));

    params.delta = alloy(data("delta", 0.341), params.delta, content,
       bowing("delta", 0.0));
    params.gamma1 = alloy(data("gamma1", 6.98), params.gamma1, content,
       bowing("gamma1", 0.0));
    params.gamma2 = alloy(data("gamma2", 2.06), params.gamma2, content,
       bowing("gamma2", 0.0));
    params.gamma3 = alloy(data("gamma3", 2.93), params.gamma3, content,
       bowing("gamma3", 0.0));

    params.def_vol_X = alloy(data("abs_def_pot_X", -0.16),
        params.def_vol_X, content, bowing("abs_def_pot_X", 0.0));
    params.def_uniax_X = alloy(data("uniax_def_pot_X", 14.26),
        params.def_uniax_X, content, bowing("uniax_def_pot_X", 0.0));
    params.def_vol_L = alloy(data("abs_def_pot_L", -4.91),
        params.def_vol_L, content, bowing("abs_def_pot_L", 0.0));
    params.def_uniax_L = alloy(data("uniax_def_pot_L", 6.5),
        params.def_uniax_L, content, bowing("uniax_def_pot_L", 0.0));

    permittivity = alloy(data("permittivity", 12.93), permittivity, content,
       bowing("permittivity", 0.0));
    _e_mobility = alloy(data("electron_mobility", 1000.0),
        _e_mobility, content, bowing("electron_mobility", 0.0));
    _h_mobility = alloy(data("hole_mobility", 200.0),
        _h_mobility, content, bowing("hole_mobility", 0.0));
  }
  else
  {
    WzDDsemiconductor* sc = dynamic_cast<WzDDsemiconductor*>(_bulk_model);
    WzDDsemiconductor::WzDDparameters& params = sc->get_parameters();

    params.EgGamma = alloy(data("Eg_G", 3.51), params.EgGamma, content,
        bowing("Eg_G", 0.0));
    params.Ev = alloy(data("E_v", -0.726), params.Ev, content,
        bowing("E_v", 0.0));

    params.m_c_zz = alloy(data("m_c_zz", 0.20), params.m_c_zz, content,
        bowing("m_c_zz", 0.0));
    params.m_c_xx = alloy(data("m_c_xx", 0.20), params.m_c_xx, content,
        bowing("m_c_xx", 0.0));

    params.A1 = alloy(data("A1", -7.21), params.A1, content,
        bowing("A1", 0.0));
    params.A2 = alloy(data("A2", -0.44), params.A2, content,
        bowing("A2", 0.0));
    params.A3 = alloy(data("A3", 6.68), params.A3, content,
        bowing("A3", 0.0));
    params.A4 = alloy(data("A4", -3.46), params.A4,content, 
        bowing("A4", 0.0));
    params.A5 = alloy(data("A5", -3.40), params.A5, content,
        bowing("A5", 0.0));
    params.A6 = alloy(data("A6", -4.90), params.A6,content, 
        bowing("A6", 0.0));

    params.a_x = alloy(data("a_x", -4.9), params.a_x, content,
        bowing("a_x", 0.0));
    params.a_z = alloy(data("a_z", -11.3), params.a_z,content, 
        bowing("a_z", 0.0));

    params.D1 = alloy(data("D1", -3.7), params.D1, content,
        bowing("D1", 0.0));
    params.D2 = alloy(data("D2", 4.5), params.D2,content, 
        bowing("D2", 0.0));
    params.D3 = alloy(data("D3", 8.2), params.D3, content,
        bowing("D3", 0.0));
    params.D4 = alloy(data("D4", -4.1), params.D4, content,
        bowing("D4", 0.0));
    params.D5 = alloy(data("D5", -4.0), params.D5, content,
        bowing("D5", 0.0));
    params.D6 = alloy(data("D6", -5.5), params.D6, content,
        bowing("D6", 0.0));
    params.delta_s = alloy(data("delta_s", 0.017), params.delta_s,
        content, bowing("delta_s", 0.0));
    params.delta_cr = alloy(data("delta_cr", 0.010), params.delta_cr,
        content, bowing("delta_cr", 0.0));

    permittivity = alloy(data("permittivity", 9.5), permittivity,
        content, bowing("permittivity", 0.0));
    _e_mobility = alloy(data("electron_mobility", 1000.0),
      _e_mobility, content, bowing("electron_mobility", 0.0));
    _h_mobility = alloy(data("hole_mobility", 200.0),
      _h_mobility, content, bowing("hole_mobility", 0.0));
  }
}

void
SemiconductorModel::prepare_element_data(void)
{
  if (!_is_prepared)
  {
    try
    {
      calculate_equilibrium_properties(BOTH, SimulationOptions::T);
    }
    catch (...)
    {
    }
  
    _is_prepared = true;
  }
}


//
// TODO
// Very crude implementation at the moment
// 
void
SemiconductorModel::extract_band_properties(void)
{
  // treat conduction band
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    _bulk_model->get_conduction_band_energy_mass();
  // get minimum
  int id = 0;
  for (int i = 1; i < cbs.size(); i++)
  {
    if (cbs[i].energy < cbs[id].energy)
      id = i;
  }
  get_conduction_band().band_edge = cbs[id].energy;
  get_conduction_band().effective_mass = cbs[id].mass_DOS
    * std::pow(cbs[id].degeneracy, 2.0 / 3.0);
  
  // treat valence band
  const std::vector<DDsemiconductor::band_extremum>& vbs =
    _bulk_model->get_valence_band_energy_mass();
  // get maximum
  id = 0;
  double kT = SimulationOptions::T * Constants::k_B;
  double delta_max = 4.0 * kT;
  for (int i = 1; i < vbs.size(); i++)
  {
    if (vbs[i].energy > vbs[id].energy)
      id = i;
  }
  get_valence_band().band_edge = vbs[id].energy;
  double tmp = 0;
  // include other bands
  for (int i = 0; i < vbs.size(); i++)
  {
    double delta = get_valence_band().band_edge - vbs[i].energy;
    if (delta < delta_max)
      tmp += vbs[i].degeneracy * std::pow(vbs[i].mass_DOS, 1.5)
        * std::exp(-delta / kT);
  }
  get_valence_band().effective_mass = std::pow(tmp, 2.0 / 3.0);
  
}

void
SemiconductorModel::print_info(void) const
{
  _bulk_model->calculate_conduction_band_extremum();
  const std::vector<DDsemiconductor::band_extremum>& cbs =
    _bulk_model->get_conduction_band_energy_mass();
  cout << " - conduction bands:\n";
  for (int i = 0 ; i < cbs.size(); i++)
  {
    cout << "   Ec = " << cbs[i].energy
      << ", m = " << cbs[i].mass_DOS
      << ", d = " << cbs[i].degeneracy << endl;
  }
  cout << "   Nc = " << get_conduction_band().effective_DOS << " cm^-3"
    << "  m_dos = " << get_conduction_band().effective_mass << "\n";

  _bulk_model->calculate_valence_band_extremum();
  const std::vector<DDsemiconductor::band_extremum>& vbs =
    _bulk_model->get_valence_band_energy_mass();
  cout << " - valence bands:\n";
  for (int i = 0 ; i < vbs.size(); i++)
  {
    cout << "   Ev = " << vbs[i].energy
      << ", m = " << vbs[i].mass_DOS
      << ", d = " << vbs[i].degeneracy << endl;
  }
  cout << "   Nv = " << get_valence_band().effective_DOS << " cm^-3"
    << "  m_dos = " << get_valence_band().effective_mass << "\n";

  cout << " - Ef0 = " << get_equilibrium_fermi_level()
    << ", ni^2 = " << get_intrinsic_density_squared();
  cout << endl;
}

void
SemiconductorModel::calculate_equilibrium_properties(int coupling,
    double temperature)
{

  assert(_bulk_model != NULL);

  // calculate conduction and valence band data
  _bulk_model->calculate_conduction_band_extremum();
  _bulk_model->calculate_valence_band_extremum();

  // get the band properties from _bulk_model
  extract_band_properties();

  // call the method of the parent class
  Parent::calculate_equilibrium_properties(coupling, temperature);

}


void
SemiconductorModel::calculate_all(double potential,
    double fermi_e, double fermi_h, const Point& coord)
{
  int coupling = get_coupling_type();

  // in this simple model all temperatures are equal
  double kT = electron_vt;

  // call the method of the parent class
  Parent::calculate_all(potential, fermi_e, fermi_h, coord);

  double n = electron_density;
  double dn = electron_density_derivative;
  double p = hole_density;
  double dp = hole_density_derivative;

  // 4.) mobilities / conductivities
  // For both statistics:
  // 
  //   mu_n = e * D_n * (1 / n) * (dn / dEf_e)
  //
  // NOTE: kT := kB * T / e includes already e
  //if (coupling & DriftDiffusionDefs::ELECTRONS)
  //{
    electron_mobility = _e_mobility;
    electron_conductivity = _e_mobility * n;
    electron_conductivity_derivatives[0] = _e_mobility * dn;
    electron_conductivity_derivatives[1] = _e_mobility * dn;
  //}
  //if (coupling & DriftDiffusionDefs::HOLES)
  //{
    hole_mobility = _h_mobility;
    hole_conductivity = _h_mobility * p;
    hole_conductivity_derivatives[0] = _h_mobility * dp;
    hole_conductivity_derivatives[2] = _h_mobility * dp;
  //}

    if (_exciton_gen_param > 1e-56)
      calculate_exciton_generation();
}

