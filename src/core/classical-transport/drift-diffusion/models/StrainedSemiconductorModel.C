// $Id$

#include "StrainedSemiconductorModel.h"
#include "DDsemiconductor.h"
#include "Macrostrain.h"
#include "Utils.h"

#include "elem.h"

#include <iostream>
#include <map>
#include <vector>


TIBER_MODULE(StrainedSemiconductorModel, strained)



using namespace std;
using namespace DriftDiffusionDefs;

StrainedSemiconductorModel::StrainedSemiconductorModel(void)
  : strain_model_(NULL),
    ignore_strain_(false)
{
}



void
StrainedSemiconductorModel::prepare_element_data(void)
{
  if (!ignore_strain_)
  {
    const Elem* elem = get_element();
    assert(elem != NULL);

    const DataMap::const_iterator end = get_data_map().end();
    const DataMap::const_iterator it = get_data_map().find(elem);
    if ((it == end) || _recompute_band_parameters)
    {
      // where to put the elemental data
      ElementData& elem_data = get_data_map()[elem];
      
      /*
      map<ID, double> data;
      bool ok = strain_model_->get_solution(elem, elem->centroid(),
          _strain_ids_set, data);

      if (ok)
      {
        // set strain
        get_strain()(1,1) = data[_strain_ids[0]];
        get_strain()(2,2) = data[_strain_ids[1]];
        get_strain()(3,3) = data[_strain_ids[2]];
        get_strain()(2,1) = data[_strain_ids[3]];
        get_strain()(3,2) = data[_strain_ids[4]];
        get_strain()(3,1) = data[_strain_ids[5]];

        elem_data.polarization(0) = data[_strain_ids[6]]; 
        elem_data.polarization(1) = data[_strain_ids[7]]; 
        elem_data.polarization(2) = data[_strain_ids[8]]; 
        set_polarization(elem_data.polarization);
      }
      else
        get_strain() = 0;
      */

      get_physical_model()->set_strain(get_strain());

      // call method of parent class
      SemiconductorModel::calculate_equilibrium_properties();

      // put them into get_data_map()
      elem_data.Ec = get_conduction_band_edge();
      elem_data.Ev = get_valence_band_edge();
      elem_data.mc = get_conduction_band().effective_mass;
      elem_data.mv = get_valence_band().effective_mass;
      elem_data.Ef0 = get_equilibrium_fermi_level();
      elem_data.ni = get_intrinsic_density();

    }
    else
    {
      const ElementData& elem_data = it->second;

      set_polarization(elem_data.polarization);

      get_conduction_band().band_edge = elem_data.Ec; 
      get_conduction_band().effective_mass = elem_data.mc; 
      get_valence_band().band_edge = elem_data.Ev; 
      get_valence_band().effective_mass = elem_data.mv; 

      equilibrium_fermi_level = elem_data.Ef0;
      intrinsic_density = elem_data.ni;

      // this sets the band edges and the effective DOS in the base class
      setup_band_edges();
    }

  }
  else
    SemiconductorModel::prepare_element_data();
}


void
StrainedSemiconductorModel::copy_from(const PhysicalModelInterface* rhs)
{
  SemiconductorModel::copy_from(rhs);

  const StrainedSemiconductorModel* sc =
    dynamic_cast<const StrainedSemiconductorModel*>(rhs);

  strain_model_ = sc->strain_model_;
  ignore_strain_ = sc->ignore_strain_;
}


void
StrainedSemiconductorModel::do_init(void)
{
  SemiconductorModel::do_init();

  
  string strain_sim = get_parameter("strain_simulation", "");

  if (strain_sim == "")
    ignore_strain_ = true;
  else
  {
    // find the strain calculation to use
    strain_model_ = SimulationInterface::find_simulation(strain_sim);

    if (strain_model_ == NULL)
    {
      string msg("Simulation " + string(strain_sim) + " not found,");
      msg += " but needed for StrainedSemiconductorModel.";
      throw InitFailedException(msg);
    }

    _strain_ids.resize(9);
    _strain_ids[0] = strain_model_->get_variable_id("eps_xx");
    _strain_ids[1] = strain_model_->get_variable_id("eps_yy");
    _strain_ids[2] = strain_model_->get_variable_id("eps_zz");
    _strain_ids[3] = strain_model_->get_variable_id("eps_xy");
    _strain_ids[4] = strain_model_->get_variable_id("eps_yz");
    _strain_ids[5] = strain_model_->get_variable_id("eps_xz");
    _strain_ids[6] = strain_model_->get_variable_id("Px");
    _strain_ids[7] = strain_model_->get_variable_id("Py");
    _strain_ids[8] = strain_model_->get_variable_id("Pz");

    for (unsigned int i = 0; i < 9; i++)
    {
      // check
      if (_strain_ids[i] == INVALID_ID)
      {
        string msg("Simulation " + string(strain_sim));
        msg += " does not provide all variables needed";
        msg += " for StrainedSemiconductorModel.";
        throw InitFailedException(msg);
      }
      _strain_ids_set.insert(_strain_ids[i]);
    }

    
  }

}


void
StrainedSemiconductorModel::do_print_info(void)
{
  string space("    ");

  cout << space << "strained semiconductor model";
  if (ignore_strain_)
    cout << " (but strain is ignored!!)";
  cout << endl;
  if (!ignore_strain_ && (strain_model_ != NULL))
    cout << space << "strain simulation: " << strain_model_->get_name() << endl;
  if (_recompute_band_parameters)
    cout << "recompute band parameters before each simulation" << endl;
  
  if (SimulationOptions::verbose() > 1)
  {
    cout << endl;
    cout << space << "unstrained band parameters" << endl;
    set_lattice_temperature(SimulationOptions::T);
    get_physical_model()->calculate_conduction_band_extremum();
    get_physical_model()->calculate_valence_band_extremum();
    extract_band_properties();
    setup_band_edges();
    calculate_equilibrium_properties();

    double deg = std::pow(2.0, 2.0 / 3.0);

    const std::vector<DDsemiconductor::band_extremum>& cbs =
      get_physical_model()->get_conduction_band_energy_mass();
    cout << space << " - conduction bands:\n";
    for (int i = 0 ; i < cbs.size(); i++)
    {
      cout << space << "   Ec = " << cbs[i].energy
        << ", m = " << cbs[i].mass_DOS
        << ", d = " << cbs[i].degeneracy << endl;
    }
    cout << space << "   Nc = " << get_conduction_band().effective_DOS << " cm^-3"
      << "  m_dos = " << get_conduction_band().effective_mass / deg << "\n";

    //_bulk_model->calculate_valence_band_extremum();
    const std::vector<DDsemiconductor::band_extremum>& vbs =
      get_physical_model()->get_valence_band_energy_mass();
    cout << space << " - valence bands:\n";
    for (int i = 0 ; i < vbs.size(); i++)
    {
      cout << space << "   Ev = " << vbs[i].energy
        << ", m = " << vbs[i].mass_DOS
        << ", d = " << vbs[i].degeneracy << endl;
    }
    cout << space << "   Nv = " << get_valence_band().effective_DOS << " cm^-3"
      << "  m_dos = " << get_valence_band().effective_mass / deg << "\n";

    cout << space << " - Eg = " <<
      get_conduction_band().band_edge - get_valence_band().band_edge <<
      ", Ef0 = " << get_equilibrium_fermi_level()
      << ", ni = " << std::sqrt(get_intrinsic_density_squared());
    cout << endl;
  }
}
