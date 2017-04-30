// $Id: InterpolatedMobility.C 4135 2015-09-25 10:19:38Z maufder $

#include "InterpolatedMobility.h"
#include "DriftDiffusionProperties.h"
#include "SimulationInterface.h"
#include "Database.h"
#include "Messages.h"
#include "Constants.h"

#include "TiberModule.h"

using namespace std;


void
InterpolatedMobility::read_database(void)
{
  const Database& db = get_database();
}



void
InterpolatedMobility::do_init(void)
{
  string sim = get_option("interpolation_module", "");
  _interpolation_sim = SimulationInterface::find_simulation(sim);

  if (_interpolation_sim == nullptr)
  {
    string msg("Interpolation module '" + string(sim) + "' not found");
    throw InitFailedException(msg);
  }

  _model_name = get_option("model_name", "");
  _temperature_var = get_option("temperature_var", "");
  _field_var = get_option("field_var", "");
  _density_var = get_option("density_var", "");

  _model_id = _interpolation_sim->get_value_id(_model_name);
  _temperature_id = _interpolation_sim->get_param_id(_temperature_var, _model_id);
  _field_id = _interpolation_sim->get_param_id(_field_var, _model_id);
  _density_id = _interpolation_sim->get_param_id(_density_var, _model_id);
}



double
InterpolatedMobility::get_mobility(void)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature() / Constants::kb;
  const libMesh::RealGradient& F = dd.get_electric_field();
  double n  = dd.get_q_density(get_carrier());

  std::map<ID, double> params;
  params.insert( make_pair(_temperature_id, T) );
  params.insert( make_pair(_field_id, F.size()) );
  params.insert( make_pair(_density_id, n) );

  double mu = _interpolation_sim->get_value_and_derivative(_model_id, params).first;
  return mu;
}


double
InterpolatedMobility::get_derivative_potential(void)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature() / Constants::kb;
  const libMesh::RealGradient& F = dd.get_electric_field();
  double n  = dd.get_q_density(get_carrier());
  double dn = dd.get_q_density_derivative(get_carrier());

 // if (n < 1e15)
   // return 0.0;

  std::map<ID, double> params;
  params.insert( make_pair(_temperature_id, T) );
  params.insert( make_pair(_field_id, F.size()) );
  params.insert( make_pair(_density_id, n) );

  double dmu_dn = _interpolation_sim->get_value_and_derivative(_model_id, params, _density_id).second;

  //if (!isfinite(dmu_dn))
    //dmu_dn = 0.0;

  return dmu_dn * dn;
}

void
InterpolatedMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}

void
InterpolatedMobility::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double T = dd.get_lattice_temperature() / Constants::kb;
  const libMesh::RealGradient& F = dd.get_electric_field();
  double n  = dd.get_q_density(get_carrier());

  std::map<ID, double> params;
  params.insert( make_pair(_temperature_id, T) );
  params.insert( make_pair(_field_id, F.size()) );
  params.insert( make_pair(_density_id, n) );


  double dmu_dF = _interpolation_sim->get_value_and_derivative(_model_id, params, _field_id).second;
  double fac = (F.size() > 0) ? dmu_dF / F.size() : 0.0;
  dm = -F * fac;

  //if (!isfinite(dmu_dF) || !isfinite(fac))
    //dm.zero();

 // if (n < 1e15)
   // dm.zero();
}

void
InterpolatedMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
}

