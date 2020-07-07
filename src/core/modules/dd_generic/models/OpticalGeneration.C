// $Id: OpticalGeneration.C 4135 2015-09-25 10:19:38Z maufder $


#include "OpticalGeneration.h"
#include "Material.h"
#include "SimulationInterface.h"
#include "DriftDiffusionProperties.h"
#include "TiberMath.h"

#include "TiberModule.h"

#include <string>
#include "elem.h"

// add
#include "ExtProfile1D.h"
#include "ExternalProfile.h"


using namespace std;

void
OpticalGeneration::do_init(void)
{
  RecombinationModelInterface::do_init();

  string gen_str(get_option("generation", "0"));
  istringstream is(gen_str);

  if ((get_carrier_names().size() != 2) && (get_carrier_names().size() != 1))
    throw InitFailedException("Optical generation model needs exactly "
        "one or two carriers");

  get_parameter("multiplier", _multiplier);
  get_parameter("use_occupation_factors", _use_occupation);
 
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double val;

  if (get_options().has_submodel("profile"))
  {
    _profile = ExternalProfile::create(get_options().submodels_begin("profile")->second);
  }
  else if ((is >> val) || (gen_str[0] == '$'))
  {
    get_parameter("generation", _generation);
  }
  else 
  {
    vector<string> gens;
    Utils::extract_vector(gen_str, gens);
    _generation_model.resize(gens.size());
    _gen_id.resize(gens.size());
    for (size_t i = 0; i < gens.size(); ++i)
    {
      _generation_model[i] = SimulationInterface::find_simulation(gens[i]);
      if (_generation_model[i] == NULL)
        throw InitFailedException("Cannot find generation model: " + gens[i]);

      _gen_id[i] = _generation_model[i]->get_solution_id("Generation");
    }
  }
}

void
OpticalGeneration::calculate_rate_and_derivatives(std::vector<double>& R,
    std::vector<std::vector<double>>& dPotentials)
{

  if (_generation_model.size() > 0)
  {

    DriftDiffusionProperties& dd = get_driftdiffusionproperties();

    const Elem* el = dd.get_element();

    _generation = 0;

    vector<double> tmp(1);
    for (size_t i = 0; i < _generation_model.size(); ++i)
    {
      if (_generation_model[i]->get_solution(el, _gen_id[i], tmp,
          vector<Point>(1, dd.get_coordinates())))
        _generation += tmp[0];
    }

  }
  else if  (_profile  != nullptr)
  {
    DriftDiffusionProperties& dd = get_driftdiffusionproperties();
    _generation = _profile->get_data(dd.get_element(), dd.get_coordinates());
  }

  
  double rate = _multiplier * _generation;


  if (get_carrier_names().size() == 1)
  {
    ID id1 = this->get_carrier_ids()[0];
    const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

    double kT = dd.get_lattice_temperature();

    double E1  = dd.get_carrier_properties(id1)->get_band_edge();

    double ch1 = dd.get_carrier_properties(id1)->get_charge();

    double Ef1 = -dd.get_q_fermi_potential(id1);

    double f1, df1;
    if (ch1 > 0)
    {
      auto ff = Distributions::fermi_dirac(-Ef1 + E1, kT);
      f1 = ff.first;
      df1 = -ff.second;
    }
    else
    {
      auto ff = Distributions::fermi_dirac(Ef1 - E1, kT);
      f1 = ff.first;
      df1 = ff.second;
    }

    double g1 = 1 - f1;
    double dg1 = -df1;

    //R[id1] = -(rate * g1);
    R[id1] = -rate;
    dPotentials[id1][id1] = 0;
    dPotentials[id1][dd.n_known_carriers()] = 0;
    //dPotentials[id1][id1] = (rate * dg1);
    //dPotentials[id1][dd.n_known_carriers()] = ch1 * (rate * dg1);
  }
  else
  {
    // energies are referred to negatively charged particles, so
    // taking away one negative carrier means adding one 'hole'

    const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

    ID id1 = this->get_carrier_ids()[0];
    ID id2 = this->get_carrier_ids()[1];


    double kT = dd.get_lattice_temperature();

    double E1  = dd.get_carrier_properties(id1)->get_band_edge();
    double E2  = dd.get_carrier_properties(id2)->get_band_edge();

    if (E1 > E2)
    {
      swap(id1, id2);
      swap(E1, E2);
    }

    double ch1 = dd.get_carrier_properties(id1)->get_charge();
    double ch2 = dd.get_carrier_properties(id2)->get_charge();

    double n1  = dd.get_q_density(id1);
    double n2  = dd.get_q_density(id2);
    double N1  = dd.get_carrier_properties(id1)->get_maximum_density();
    double N2  = dd.get_carrier_properties(id2)->get_maximum_density();
    double dn1 = dd.get_q_density_derivative(id1);
    double dn2 = dd.get_q_density_derivative(id2);


    double g1 = 1, g2 = 1, dg1 = 0, dg2 = 0;
    double sign1 = 1, sign2 = 1;

    if (ch1 > 0)
    {
      // it's a hole
      g1 = 1 - n1/N1;
      dg1 = -dn1/N1;
      sign1 = -1;
    }
    else
    {
      g1 = n1/N1;
      dg1 = dn1/N1;
    }

    if (ch2 > 0)
    {
      g2 = n2/N2;
      dg2 = dn2/N2;
    }
    else
    {
      g2 = 1 - n2/N2;
      dg2 = -dn2/N2;
      sign2 = -1;
    }

    if (!_use_occupation)
    {
      g1 = g2 = 1;
      dg1 = dg2 = 0;
    }

    R[id1] = sign1 * rate * g1 * g2;
    R[id2] = sign2 * rate * g1 * g2;
    //cerr << "   " << R[id1] << " " << R[id2] << endl;

    double dR1 = rate * dg1 * g2;
    double dR2 = rate * g1 * dg2;
    dPotentials[id1][id1] = -sign1 * dR1;
    dPotentials[id1][id2] = -sign1 * dR2;
    dPotentials[id2][id1] = -sign2 * dR1;
    dPotentials[id2][id2] = -sign2 * dR2;
    dPotentials[id1][dd.n_known_carriers()] = sign1 * (dR1 + dR2);
    dPotentials[id2][dd.n_known_carriers()] = sign2 * (dR1 + dR2);
  }

}



