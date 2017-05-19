
#include "ExcitonEnergyTransfer.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"


using namespace std;

ExcitonEnergyTransfer::QRecMap
ExcitonEnergyTransfer::_qrec_vals;

void
ExcitonEnergyTransfer::read_database(void)
{
  const Database& db = get_database();

  db.set_section("permittivity");
  _er = db.get("permittivity", 1.0);
}



void
ExcitonEnergyTransfer::do_init(void)
{
  RecombinationModelInterface::do_init();

  _tau = get_option("tau", _tau);
  _tau_rad = get_option("tau", _tau_rad);
  _Rf   = get_option("forster_radius", _Rf);
  _Rd   = get_option("dexter_radius", _Rd);
  _R_da = get_option("average_distance", _R_da);

  if (get_carrier_names().size() != 2)
    throw InitFailedException("Exciton EnergyTransfer model needs exactly two recombining carriers");

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  for (auto name : get_carrier_names())
  {
    if (!dd.get_carrier_properties(name)->is_exciton())
      throw InitFailedException("Recombination '" + get_default_name() + ": carrier '" + name + "' is not an exciton");
  }

  _id_d = this->get_carrier_ids()[0];
  _id_a = this->get_carrier_ids()[1];

  double E1  = dd.get_carrier_properties(_id_d)->get_band_edge();
  double E2  = dd.get_carrier_properties(_id_a)->get_band_edge();

  if (E1<E2)
    swap(_id_d, _id_a);

  _m = dd.get_carrier_properties(_id_d)->get_effective_mass();

  double RBeff = Constants::bohr_radius * 100.0 * _er / _m;  //effective Bohr radius in cm

  _Kf = (get_option("forster", true)) ? pow(_Rf / _R_da, 6.0) / _tau_rad            : 0.0; //Forster rate
  _Kd = (get_option("dexter", true) ) ? 1e3 * exp(2.0*(_Rd - _R_da) / RBeff) / _tau : 0.0; //Dexter rate

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
ExcitonEnergyTransfer::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double Ed = dd.get_carrier_properties(_id_d)->get_band_edge();
  double Ea = dd.get_carrier_properties(_id_a)->get_band_edge();

  double kT = dd.get_lattice_temperature();
  double xd  = dd.get_q_density(_id_d);
  double xa  = dd.get_q_density(_id_a);
  double Na  = dd.get_carrier_properties(_id_a)->get_effective_DOS();
  double dxd = dd.get_q_density_derivative(_id_d);
  double dxa = dd.get_q_density_derivative(_id_a);
  double fd = dd.get_q_fermi_potential(_id_d);
  double fa = dd.get_q_fermi_potential(_id_a);

  double beta = 1.0/kT;
  double exponential = exp( beta*(fd-fa) );
  double stat = 1.0 - exponential; 

  double rate = (_Kf + _Kd) * stat * xd * (Na + xa) / Na;
  double derd = - (_Kf + _Kd) * (Na + xa) * (dxd * stat + beta * xd * exponential) / Na;
  double dera = - (_Kf + _Kd) * xd * (dxa * stat + beta * (Na + xa) * exponential) / Na;

  R[_id_d] = rate;
  R[_id_a] = -rate;

  dPotentials[_id_d][_id_d] =  derd;
  dPotentials[_id_d][_id_a] =  dera;
  dPotentials[_id_a][_id_d] = -derd;
  dPotentials[_id_a][_id_a] = -dera;

}

void
ExcitonEnergyTransfer::do_reinit(void)
{

}
