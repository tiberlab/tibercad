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
        "two recombining carriers");

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

  get_parameter("extraction_barrier", _extraction_barrier);

}




void
DirectRecombination::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{

  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const char ct1 = dd.get_carrier_properties(id1)->get_carrier_type();
  const char ct2 = dd.get_carrier_properties(id2)->get_carrier_type();

  double kT = dd.get_lattice_temperature();

  if (ct1 != ct2)
  {
    // this is the standard direct recombination for electrons hole pairs
    if (ct1 != 'e')
      swap(id1, id2);

    double Ef1 = -dd.get_q_fermi_potential(id1);
    double Ef2 = -dd.get_q_fermi_potential(id2);
    double beta = 1.0/kT;

    double n1  = dd.get_q_density(id1);
    double n2  = dd.get_q_density(id2);
    double dn1  = dd.get_q_density_derivative(id1);
    double dn2  = dd.get_q_density_derivative(id2);
    double q1 = dd.get_carrier_properties(id1)->get_charge();
    double q2 = dd.get_carrier_properties(id2)->get_charge();


    //double E01 = dd.get_carrier_properties(id1)->get_band_edge();
    //double E02 = dd.get_carrier_properties(id2)->get_band_edge();

    long double arg = (Ef2 - Ef1) * beta;
    long double exponential = 0;
    long double stat_fac = 1.0;
    if (fabs(arg) < 1e-3)
    {
      exponential = 1 + arg;
      stat_fac = -arg;
    }
    else
    {
      exponential = exp(arg);
      stat_fac = 1.0 - exponential;
    }

    long double g = C_ * n1 * n2;

    R[id1] = g * stat_fac;
    R[id2] = g * stat_fac;
    //cerr << "n1 = " << n1 << " n2 = " << n2 << " -> " << R[id1] << endl;

    long double dR0 = stat_fac * C_ * (n2 * dn1 + n1 * dn2);
    long double dR1 = -C_ * n2 * (dn1 * stat_fac + beta * n1 * exponential);
    long double dR2 = -C_ * n1 * (dn2 * stat_fac - beta * n2 * exponential);

    dPotentials[id1][id1] = dR1;
    dPotentials[id1][id2] = dR2;
    dPotentials[id2][id1] = dR1;
    dPotentials[id2][id2] = dR2;
    dPotentials[id1][dd.n_known_carriers()] = dR0;
    dPotentials[id2][dd.n_known_carriers()] = dR0;
  }
  else
  {
    double E1  = dd.get_carrier_properties(id1)->get_band_edge();
    double E2  = dd.get_carrier_properties(id2)->get_band_edge();
    double q1  = dd.get_carrier_properties(id1)->get_charge();
    double q2  = dd.get_carrier_properties(id2)->get_charge();

    if ((ct1 == 'e') && (E2 < E1))
    {
      swap(id1, id2);
      swap(E1, E2);
    }
    else if ((ct1 == 'h') && (E1 < E2))
    {
      swap(id1, id2);
      swap(E1, E2);
    }


    double n1  = dd.get_q_density(id1);
    double n2  = dd.get_q_density(id2);
    double N1  = dd.get_carrier_properties(id1)->get_maximum_density();
    double N2  = dd.get_carrier_properties(id2)->get_maximum_density();
    double dn1 = dd.get_q_density_derivative(id1);
    double dn2 = dd.get_q_density_derivative(id2);
    double Ef1 = -dd.get_q_fermi_potential(id1);
    double Ef2 = -dd.get_q_fermi_potential(id2);

    if (_extraction_barrier > 0)
      kT = _extraction_barrier;

    double beta = (ct1 == 'e') ? 1.0/kT : -1.0/kT;

    double thermal = exp(-fabs(E2 - E1) / kT);
    double arg1 = (Ef2 - Ef1) * beta;
    double arg2 = -arg1;

    double exponential1 = 0;
    double stat1 = 1.0;
    if (fabs(arg1) > 1e-3)
    {
      exponential1 = exp(arg1);
      stat1 = 1.0 - exponential1;
    }
    else
    {
      exponential1 = 1 + arg1;
      stat1 = -arg1;
    }

    double exponential2 = 0;
    double stat2 = 1.0;
    if (fabs(arg2) > 1e-3)
    {
      exponential2 = exp(arg2);
      stat2 = 1.0 - exponential2;
    }
    else
    {
      exponential2 = 1 + arg2;
      stat2 = -arg2;
    }
    double n1f = n1 / N1;
    double n2f = n2 / N2;

    double C1 = 0.5*C_*thermal;
    double C2 = 0.5*C_;

    R[id1] = C1 * stat1 * n1 * (1.0 - n2f) - C2 * stat2 * n2 * (1.0 - n1f);
    R[id2] = -R[id1];


    double dR0 =  C1 * stat1 * ( (1 - n2f)*dn1 - n1*dn2/N2 )
        - C2 * stat2 * ( (1 - n1f)*dn2 - n2 * dn1/N1 );
    double dR1 = -C1 * (1 - n2f) * ( dn1 * stat1 + beta * n1 * exponential1)
        - C2 * n2 * ( dn1/N1 * stat2 + beta * (1 - n1f) * exponential2);
    double dR2 = -C1 * n1 * (-dn2/N2 * stat1 - beta * (1 - n2f) * exponential1)
        + C2 * (1 - n1f) * (dn2 * stat2 + beta * n2 * exponential2);
    //cerr << R[id1] << " : " << dR0 <<  " " << dR1 << " " << dR2 << endl;

    dPotentials[id1][id1] =  dR1;
    dPotentials[id1][id2] =  dR2;
    dPotentials[id2][id1] = -dR1;
    dPotentials[id2][id2] = -dR2;
    dPotentials[id1][dd.n_known_carriers()] =  dR0;
    dPotentials[id2][dd.n_known_carriers()] = -dR0;
  }
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
