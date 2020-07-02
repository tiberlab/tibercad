// $Id: AvalancheGeneration.C 4145 2015-10-02 11:53:20Z maufder $

#include "AvalancheGeneration.h"
#include "DriftDiffusionProperties.h"
#include "TiberMath.h"

#include "Material.h"
#include "Database.h"


#include "TiberModule.h"


using namespace std;

AvalancheGeneration::AvalancheGeneration(const ModelOptions& options) :
  RecombinationModelInterface(options),
  _w0(0.05)
{
}



AvalancheGeneration::~AvalancheGeneration(void)
{
}


void
AvalancheGeneration::read_database(void)
{
  const Database& db = get_database();

  db.set_section("recombination/AvalancheGeneration");

  _w0 = db.get("w0", _w0);

  db.get("a", _a_param);
  db.get("b", _b_param);
}




void
AvalancheGeneration::do_init(void)
{
  RecombinationModelInterface::do_init();

  unsigned int nc = get_carrier_names().size();

  _a_param.resize(nc, 0);
  _b_param.resize(nc, 0);

  get_option("a", _a_param);
  get_option("b", _b_param);
  get_parameter("w0", _w0);

}



void
AvalancheGeneration::calculate_rate_and_derivatives(std::vector<double>& R,
    std::vector<std::vector<double>>& dPotentials)
{
  const vector<ID>& ids = this->get_carrier_ids();
  unsigned int nc = ids.size();

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  //double Efn = -dd.get_q_fermi_potential(id1);
  //double Efp = -dd.get_q_fermi_potential(id2);
  double kT = dd.get_lattice_temperature();

  vector<double> alpha(nc, 0.0);

  // we add a small value to be able to divide by F in any case
  double F = dd.get_electric_field().norm() + 1e-6;

  for (unsigned int i = 0; i < nc; ++i)
  {
    ID id = ids[i];

    libMesh::RealGradient grad;
    dd.get_grad_fermi(id, grad);

    F = grad.norm() + 1e-6;

    double arg = -_b_param[i] / F;
    if (arg < -15)
      alpha[i] = 0;
    else
      alpha[i] = _a_param[i] * exp(arg);


    double flux = dd.get_q_conductivity(id) * grad.norm();

    for (unsigned int j = 0; j < nc; ++j)
      R[ids[j]] -= alpha[i] * flux;
  }

  //dPotentials[id1][id1] = dR1;
  //dPotentials[id1][id2] = dR2;
  //dPotentials[id2][id1] = dR1;
  //dPotentials[id2][id2] = dR2;
  //dPotentials[id1][dd.n_known_carriers()] = dR0;
  //dPotentials[id2][dd.n_known_carriers()] = dR0;
}

