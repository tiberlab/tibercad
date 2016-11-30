// $Id: DirectRecombination.C 3542 2013-03-01 09:31:59Z maufder $

#include "DirectRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"


using namespace std;

DirectRecombination::QRecMap
DirectRecombination::_qrec_vals;

void
DirectRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);
}



void
DirectRecombination::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 2)
    throw InitFailedException("Direct recombination model needs exactly "
        "two recombinaing carriers");

  get_parameter("C", C_);

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
DirectRecombination::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double d0  = dd.get_q_density(get_carrier_ids()[0]);
  double d1  = dd.get_q_density(get_carrier_ids()[1]);
  double d00 = dd.get_equilibrium_q_density(get_carrier_ids()[0]);
  double d10 = dd.get_equilibrium_q_density(get_carrier_ids()[1]);
  //double ni = dd.get_intrinsic_density();
  //double gn = dd.get_electron_gamma();
  //double gp = dd.get_hole_gamma();

  recomb_e = recomb_h = C_ * (d0 * d1 - d00 * d10);
}



void
DirectRecombination::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double d0  = dd.get_q_density(get_carrier_ids()[0]);
  double d1  = dd.get_q_density(get_carrier_ids()[1]);

  recomb_e[0] = recomb_h[0] = C_ * d1; // dR/dd0
  recomb_e[1] = recomb_h[1] = C_ * d0; // dR/dd1
  recomb_e[2] = recomb_e[3] = recomb_h[2] = recomb_h[3] = 0.0;
}

double
DirectRecombination::calculate_rate_and_derivatives(std::vector<double>& dPotentials)
{
  const ID id1 = this->get_carrier_ids()[0];
  const ID id2 = this->get_carrier_ids()[1];
  if ((id1 == DriftDiffusionProperties::unknown_carrier_id) ||
      (id2 == DriftDiffusionProperties::unknown_carrier_id))
  {
    return 0.0;
  }

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Ef1 = -dd.get_q_fermi_potential(id1);
  double Ef2 = -dd.get_q_fermi_potential(id2);
  double kT = dd.get_lattice_temperature();

  double n1  = dd.get_q_density(id1);
  double n2  = dd.get_q_density(id2);
  double dn1  = dd.get_q_density_derivative(id1);
  double dn2  = dd.get_q_density_derivative(id2);
  double q1 = dd.get_carrier_properties(id1)->get_charge();
  double q2 = dd.get_carrier_properties(id2)->get_charge();

  //double E01 = dd.get_carrier_properties(id1)->get_band_edge();
  //double E02 = dd.get_carrier_properties(id2)->get_band_edge();

  double exponential = exp((Ef2 - Ef1) / kT);
  double stat_fac = 1.0 - exponential;
  double g = C_ * n1 * n2;

  double R = g * stat_fac;

  double dR1 = -C_ * n2 * (-dn1 * stat_fac + 1/kT * n1 * exponential);
  double dR2 = -C_ * n1 * (-dn2 * stat_fac - 1/kT * n2 * exponential);

  dPotentials[id1] = dR1;
  dPotentials[id2] = dR2;
  dPotentials[dd.n_known_carriers()] = stat_fac * C_ * (n2 * dn1 + n1 * dn2);

  return R;
}

void
DirectRecombination::do_reinit(void)
{
  if (_quantum_optics != NULL)
  {
    SimulationInterface* sim = SimulationInterface::get_simulation(get_simulator_id());
    unsigned int seqnum = sim->get_solve_sequence_number();

    QRecMap::iterator it(_qrec_vals.find(make_pair(_quantum_optics, sim)));
    if ((it != _qrec_vals.end()) && ((it->second).first == seqnum))
      C_ = (it->second).second;
    else
    {
      map<ID, vector<double> > data;
      data[_rec_id];

      if (_quantum_optics->get_solution(data))
      {
        double rec = data[_rec_id][0];
        Messages::info("Recalculate radiative recombination parameter: ", false);

        // now we have to integrate the term (np - ni^2)
        // for that, we have to loop over all elements
        // TODO this has to be checked for parallel execution

        const SimulationEnvironment& env = sim->get_environment();

        data.clear();
        // TODO IntrinsicDensity is currently missing
        ID edens_id = sim->get_solution_id("eDensity");
        ID hdens_id = sim->get_solution_id("hDensity");
        //ID idens_id = sim->get_solution_id("IntrinsicDensity");
        data[edens_id];
        data[hdens_id];

        unsigned int dim = sim->get_mesh().mesh_dimension();
        libMesh::UniquePtr<libMesh::FEBase> fe(sim->build_finite_element(dim, libMesh::FEType()));
        libMesh::UniquePtr<libMesh::QBase> qrule(libMesh::QBase::build(libMeshEnums::QGAUSS, dim, libMeshEnums::FIFTH));
        fe->attach_quadrature_rule(qrule.get());

        const vector<Real>& JxW = fe->get_JxW();

        double tot_rec = 0.0;

        MeshBase::const_element_iterator it = sim->active_local_elements_begin();
        const MeshBase::const_element_iterator end = sim->active_local_elements_end();
        for ( ; it != end; ++it)
        {
          const libMesh::Elem* elem = *it;
          if (_quantum_optics->includes_region(elem->subdomain_id()))
          {
            fe->reinit(elem);

            if (sim->get_solution(elem, data, qrule->get_points(), true))
            {
              for (int n = 0; n < qrule->n_points(); ++n)
              {
                double np = data[edens_id][n] * data[hdens_id][n];
                tot_rec += JxW[n] * np;
              }
            }
          }
        }

        C_ = rec / tot_rec;
        _qrec_vals[make_pair(_quantum_optics, sim)] = make_pair(seqnum, C_);

        ostringstream os;
        os << "B = " << rec / tot_rec << endl;
        Messages::info(os.str());
        Messages::newline();
      }
    }
  }
}
