// $Id$

#include "GenericRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;

namespace {
  pair<double, double> fermi_dirac(double E)
  {
    double f = 0, deriv = 0;
    double arg = E;
    if (arg > 50)
    {
      f = exp(-arg);
      deriv = -f;
    }
    else if (arg < -50)
    {
      deriv = -exp(arg);
      f = 1 + deriv;
    }
    else
    {
      double expfac = exp(arg);
      double denom = 1.0 + expfac;
      f = 1.0 / denom;
      deriv = -expfac * f / denom;
    }

    return make_pair(f, deriv);
  }
}

/*
void
GenericRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("recombination/direct");

  C_ = db.get("C", C_);
}
*/



void
GenericRecombination::do_init(void)
{
  RecombinationModelInterface::do_init();

  if (get_carrier_names().size() != 2)
    throw InitFailedException("Generic recombination model needs exactly "
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

}




void
GenericRecombination::calculate_rate_and_derivatives(std::vector<double>& R, std::vector<std::vector<double>>& dPotentials)
{

  ID id1 = this->get_carrier_ids()[0];
  ID id2 = this->get_carrier_ids()[1];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  const char ct1 = dd.get_carrier_properties(id1)->get_carrier_type();
  const char ct2 = dd.get_carrier_properties(id2)->get_carrier_type();

  double kT = dd.get_lattice_temperature();

  if (ct1 != ct2)
  {
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

    double exponential = exp((Ef2 - Ef1) * beta);
    double stat_fac = 1.0 - exponential;
    double g = C_ * n1 * n2;

    R[id1] = g * stat_fac;
    R[id2] = g * stat_fac;

    double dR0 = stat_fac * C_ * (n2 * dn1 + n1 * dn2);
    double dR1 = -C_ * n2 * (dn1 * stat_fac + beta * n1 * exponential);
    double dR2 = -C_ * n1 * (dn2 * stat_fac - beta * n2 * exponential);

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

    if ((ct1 == 'e') && (E1 < E2))
    {
      swap(id1, id2);
      swap(E1, E2);
    }
    else if ((ct1 == 'h') && (E2 < E1))
    {
      swap(id1, id2);
      swap(E1, E2);
    }

    /*
    / *if ((ct1 == 'e') && (E1 < E2))
    {
      swap(id1, id2);
      swap(E1, E2);
    }
    else if ((ct1 == 'h') && (E2 < E1))
    {
      swap(id1, id2);
      swap(E1, E2);
    }* /

    if (E1 < E2)
    {
      swap(id1, id2);
      swap(E1, E2);
    }

    // now E1 > E2 for electrons
    //     E2 > E1 for holes


    //double n1  = dd.get_q_density(id1);
    //double n2  = dd.get_q_density(id2);
    //double N1  = dd.get_carrier_properties(id1)->get_effective_DOS();
    //double N2  = dd.get_carrier_properties(id2)->get_effective_DOS();
    //double dn1 = dd.get_q_density_derivative(id1);
    //double dn2 = dd.get_q_density_derivative(id2);
    double Ef1 = -dd.get_q_fermi_potential(id1);
    double Ef2 = -dd.get_q_fermi_potential(id2);
    double beta = 1/kT; //(ct1 == 'e') ? 1.0/kT : -1.0/kT;

    double exponential = exp((Ef2 - Ef1) * beta);
    double stat = 1.0 - exponential;

    auto occ1 = fermi_dirac((E1 - Ef1) * beta);
    double f1 = occ1.first;
    double df1 = occ1.second;

    auto occ2 = fermi_dirac((Ef2 - E2) * beta);
    double f2 = occ2.first;
    double df2 = -occ2.second;

    R[id1] = C_ * f1 * f2 * stat;
    R[id2] = -R[id1];

    double dRE1 = C_ * df1 * f2 * stat;
    double dRE2 = C_ * f1 * df2 * stat;


    double dR0 =  -(dRE1 + dRE2) * beta;
    double dR1 = -dRE1 * beta + C_ * f1 * f2 * beta * exponential;
    double dR2 = -dRE2 * beta - C_ * f1 * f2 * beta * exponential;

    dPotentials[id1][id1] = - dR1;
    dPotentials[id1][id2] = - dR2;
    dPotentials[id2][id1] = dR1;
    dPotentials[id2][id2] = dR2;
    dPotentials[id1][dd.n_known_carriers()] =  dR0;
    dPotentials[id2][dd.n_known_carriers()] = -dR0;
    */

    ///*
    double n1  = dd.get_q_density(id1);
    double n2  = dd.get_q_density(id2);
    double N1  = 100 * dd.get_carrier_properties(id1)->get_effective_DOS();
    double N2  = 100 * dd.get_carrier_properties(id2)->get_effective_DOS();
    double dn1 = dd.get_q_density_derivative(id1);
    double dn2 = dd.get_q_density_derivative(id2);
    double Ef1 = -dd.get_q_fermi_potential(id1);
    double Ef2 = -dd.get_q_fermi_potential(id2);
    double beta = (ct1 == 'e') ? 1.0/kT : -1.0/kT;

    double thermal = 1; //exp( 0.5*beta*(E1-E2) - 0.5*fabs(beta*(E1-E2)) );
    double exponential = exp((Ef2 - Ef1) * beta);
    double stat = 1.0 - exponential;

    // TODO N2 should be the maximum density for n2, but now its 100*Nc
    R[id1] = C_ * thermal * stat * n1 * (N2 - n2);
    R[id2] = -R[id1];


    double dR0 =  C_ * thermal * stat * ( (N2 - n2)*dn1 - n1*dn2 );
    double dR1 = -C_ * thermal * (N2 - n2) * ( dn1 * stat + beta * n1 * exponential);
    double dR2 = -C_ * thermal * n1 * (-dn2 * stat - beta * (N2 - n2) * exponential);

    dPotentials[id1][id1] =  dR1;
    dPotentials[id1][id2] =  dR2;
    dPotentials[id2][id1] = -dR1;
    dPotentials[id2][id2] = -dR2;
    dPotentials[id1][dd.n_known_carriers()] =  dR0;
    dPotentials[id2][dd.n_known_carriers()] = -dR0;
    //*/
  }
}

void
GenericRecombination::do_reinit(void)
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
