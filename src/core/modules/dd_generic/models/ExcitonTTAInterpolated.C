
#include "ExcitonTTAInterpolated.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"


using namespace std;

ExcitonTTAInterpolated::QRecMap
ExcitonTTAInterpolated::_qrec_vals;

void
ExcitonTTAInterpolated::read_database(void)
{
  const Database& db = get_database();
}



void
ExcitonTTAInterpolated::do_init(void)
{
  RecombinationModelInterface::do_init();

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  for (auto name : get_carrier_names())
  {
    if (!dd.get_carrier_properties(name)->is_exciton())
      throw InitFailedException("Recombination '" + get_default_name() + ": carrier '" + name + "' is not an exciton");

    if (dd.get_carrier_properties(name)->get_spin() != 1.0)
      throw InitFailedException("Recombination '" + get_default_name() + ": exciton '" + name + "' is not an triplet");

  }

  string sim = get_option("interpolation_module", "");
  _interpolation_sim = SimulationInterface::find_simulation(sim);

  if (_interpolation_sim == nullptr)
  {
    string msg("Interpolation module '" + string(sim) + "' not found");
    throw InitFailedException(msg);
  }

  _model_name = get_option("model_name", "");
  _temperature_var = get_option("temperature_var", "");
  _density_var = get_option("density_var", "");

  _model_id = _interpolation_sim->get_value_id(_model_name);
  _temperature_id = _interpolation_sim->get_param_id(_temperature_var, _model_id);
  _density_id = _interpolation_sim->get_param_id(_density_var, _model_id);

  if (_model_id == INVALID_ID)
    throw InitFailedException("Model name '" + _model_name + "' not found in module '" + sim + "'");

  if (_temperature_id == INVALID_ID)
    throw InitFailedException("Variable '" + _temperature_var + "' not found in model '" + _model_name + "'");

  if (_density_id == INVALID_ID)
    throw InitFailedException("Variable '" + _density_var + "' not found in model '" + _model_name + "'");


  string quantumsim = get_option("optics_simulation", "");
  if (!quantumsim.empty())
  {
    _quantum_optics = SimulationInterface::find_simulation(quantumsim);
    if (_quantum_optics == NULL)
      throw InitFailedException("Cannot find optics simulation \'" + quantumsim + "\'");

    _rec_id = _quantum_optics->get_solution_id("Recombination");
    if (_rec_id == INVALID_ID)
      throw InitFailedException("Simulation \'" + quantumsim + "\'" +
          " does not have the needed solution \'Recombination\'");

    // TODO should check for consistency of regions
  }

}




void
ExcitonTTAInterpolated::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();


  double kT = dd.get_lattice_temperature();
  double T = kT / Constants::kb;


  for (auto id : get_carrier_ids())
  {

    double x = dd.get_q_density(id);
    double dx = dd.get_q_density_derivative(id);
    double fx = dd.get_q_fermi_potential(id);

    std::map<ID, double> params;
    params.insert( make_pair(_temperature_id, T) );
    params.insert( make_pair(_density_id, x) );

    pair<double, double> K_dK = _interpolation_sim->get_value_and_derivative(_model_id, params, _density_id);

    double beta = 1.0 /kT;
    double exponential = exp(beta * fx);
    double stat = 1.0 - exponential;

    double rate = stat * x * x * K_dK.first;

    R[id] = rate;

    dPotentials[id][id] = - x * K_dK.first * (beta * x * exponential + 2.0 * dx * stat) - x * x * stat * dx * K_dK.second;
  }

}

void
ExcitonTTAInterpolated::do_reinit(void)
{

}
