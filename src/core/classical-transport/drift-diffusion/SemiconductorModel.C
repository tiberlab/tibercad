// $Id$

#include "Alloy.h"

#include "SemiconductorModel.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "Database.h"

#include "getpot.h"

#include <iostream>

using namespace DriftDiffusionDefs;
using namespace std;


SemiconductorModel::~SemiconductorModel(void)
{
  delete _bulk_model;
}

SemiconductorModel::SemiconductorModel(void)
  : _bulk_model(NULL),
    _is_prepared(false)
{
}


void
SemiconductorModel::do_init(void)
{
  Parent::do_init();
}

/*
void
SemiconductorModel::read_bowing_parameters(void)
{
}
*/

void
SemiconductorModel::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  if (get_material()->get_structure() == "zb")
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

  
  }
}


void
SemiconductorModel::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  DriftDiffusionProperties::calculate_VCA(comp_A, comp_B, xa);

  double (*alloy)(double, double, double, double) =
    Alloy::calculate_VCA_parameter;

  const SemiconductorModel* scA =
    dynamic_cast<const SemiconductorModel*>(comp_A);
  const SemiconductorModel* scB =
    dynamic_cast<const SemiconductorModel*>(comp_B);

  if (get_material()->get_structure() == "zb")
  {
    ZbDDsemiconductor* sc = dynamic_cast<ZbDDsemiconductor*>(_bulk_model);
    ZbDDsemiconductor* sca = dynamic_cast<ZbDDsemiconductor*>(scA->_bulk_model);
    ZbDDsemiconductor* scb = dynamic_cast<ZbDDsemiconductor*>(scB->_bulk_model);
    ZbDDsemiconductor::ZbDDparameters& params = sc->get_parameters();
    const ZbDDsemiconductor::ZbDDparameters& paramsA = sca->get_parameters();
    const ZbDDsemiconductor::ZbDDparameters& paramsB = scb->get_parameters();

    params.EgGamma = alloy(paramsA.EgGamma, paramsB.EgGamma, xa, 0.0);
    params.EgL = alloy(paramsA.EgL, paramsB.EgL, xa, 0.0);
    params.EgX = alloy(paramsA.EgX, paramsB.EgX, xa, 0.0);
    params.Ev = alloy(paramsA.Ev, paramsB.Ev, xa, 0.0);

    params.m_G = alloy(paramsA.m_G, paramsB.m_G, xa, 0.0);
    params.m_t_L = alloy(paramsA.m_t_L, paramsB.m_t_L, xa, 0.0);
    params.m_l_L = alloy(paramsA.m_l_L, paramsB.m_l_L, xa, 0.0);
    params.m_t_X = alloy(paramsA.m_t_X, paramsB.m_t_X, xa, 0.0);
    params.m_l_X = alloy(paramsA.m_l_X, paramsB.m_l_X, xa, 0.0);

    params.a_c = alloy(paramsA.a_c, paramsB.a_c, xa, 0.0);
    params.a_v = alloy(paramsA.a_v, paramsB.a_v, xa, 0.0);
    params.b = alloy(paramsA.b, paramsB.b, xa, 0.0);
    params.d = alloy(paramsA.d, paramsB.d, xa, 0.0);

    params.delta = alloy(paramsA.delta, paramsB.delta, xa, 0.0);
    params.gamma1 = alloy(paramsA.gamma1, paramsB.gamma1, xa, 0.0);
    params.gamma2 = alloy(paramsA.gamma2, paramsB.gamma2, xa, 0.0);
    params.gamma3 = alloy(paramsA.gamma3, paramsB.gamma3, xa, 0.0);

    params.def_vol_X = alloy(paramsA.def_vol_X,
        paramsB.def_vol_X, xa, 0.0);
    params.def_uniax_X = alloy(paramsA.def_uniax_X,
        paramsB.def_uniax_X, xa, 0.0);
    params.def_vol_L = alloy(paramsA.def_vol_L,
        paramsB.def_vol_L, xa, 0.0);
    params.def_uniax_L = alloy(paramsA.def_uniax_L,
        paramsB.def_uniax_L, xa, 0.0);

    permittivity = alloy(scA->permittivity, scB->permittivity,
        xa, 0.0);
  }
  else
  {
    WzDDsemiconductor* sc = dynamic_cast<WzDDsemiconductor*>(_bulk_model);
    WzDDsemiconductor* sca = dynamic_cast<WzDDsemiconductor*>(scA->_bulk_model);
    WzDDsemiconductor* scb = dynamic_cast<WzDDsemiconductor*>(scB->_bulk_model);
    WzDDsemiconductor::WzDDparameters& params = sc->get_parameters();
    const WzDDsemiconductor::WzDDparameters& paramsA = sca->get_parameters();
    const WzDDsemiconductor::WzDDparameters& paramsB = scb->get_parameters();

    params.EgGamma = alloy(paramsA.EgGamma, paramsB.EgGamma, xa, 0.0);
    params.Ev = alloy(paramsA.Ev, paramsB.Ev, xa, 0.0);

    params.m_c_zz = alloy(paramsA.m_c_zz, paramsB.m_c_zz, xa, 0.0);
    params.m_c_xx = alloy(paramsA.m_c_xx, paramsB.m_c_xx, xa, 0.0);

    params.A1 = alloy(paramsA.A1, paramsB.A1, xa, 0.0);
    params.A2 = alloy(paramsA.A2, paramsB.A2, xa, 0.0);
    params.A3 = alloy(paramsA.A3, paramsB.A3, xa, 0.0);
    params.A4 = alloy(paramsA.A4, paramsB.A4, xa, 0.0);
    params.A5 = alloy(paramsA.A5, paramsB.A5, xa, 0.0);
    params.A6 = alloy(paramsA.A6, paramsB.A6, xa, 0.0); 

    params.a_x = alloy(paramsA.a_x, paramsB.a_x, xa, 0.0);
    params.a_z = alloy(paramsA.a_z, paramsB.a_z, xa, 0.0);

    params.D1 = alloy(paramsA.D1, paramsB.D1, xa, 0.0);
    params.D2 = alloy(paramsA.D2, paramsB.D2, xa, 0.0);
    params.D3 = alloy(paramsA.D3, paramsB.D3, xa, 0.0);
    params.D4 = alloy(paramsA.D4, paramsB.D4, xa, 0.0);
    params.D5 = alloy(paramsA.D5, paramsB.D5, xa, 0.0);
    params.D6 = alloy(paramsA.D6, paramsB.D6, xa, 0.0);
    params.delta_s = alloy(paramsA.delta_s, paramsB.delta_s, xa, 0.0);
    params.delta_cr = alloy(paramsA.delta_cr, paramsB.delta_cr, xa, 0.0);

    permittivity = alloy(scA->permittivity, scB->permittivity,
        xa, 0.0);

  }
}



void
SemiconductorModel::prepare_element_data(void)
{
  if (!_is_prepared)
  {
    try
    {
      calculate_equilibrium_properties();
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
SemiconductorModel::calculate_equilibrium_properties(void)
{

  assert(_bulk_model != NULL);

  // calculate conduction and valence band data
  _bulk_model->calculate_conduction_band_extremum();
  _bulk_model->calculate_valence_band_extremum();

  // get the band properties from _bulk_model
  extract_band_properties();

  // call the method of the parent class
  Parent::calculate_equilibrium_properties();

}

