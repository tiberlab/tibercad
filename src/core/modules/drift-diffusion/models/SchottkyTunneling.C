// $Id$

#include "SchottkyTunneling.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "quadrature.h"




TIBER_MODULE(SchottkyTunneling, recombination, schottky_tunneling)

using namespace std;




SchottkyTunneling::SchottkyTunneling(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _max_tunnel_length(10)
{
}



void
SchottkyTunneling::read_database(void)
{
  const Database& db = get_database();
  db.set_section("schottky_tunneling");


}



void
SchottkyTunneling::do_init(void)
{
  get_option("maximum_tunnel_length", _max_tunnel_length);
  get_option("contact_name", _contact_name);


  // calculate distance from contact and put it into
  // _elem_map if it is below the maximum tunnel length

  SimulationInterface* sim = SimulationInterface::get_simulation(get_simulator_id());
  const SimulationEnvironment& env = sim->get_environment();

  // first, build a bounding box for the contact to be faster afterwards
  BoundaryElementMap::iterator bdfirst(env.boundary_elements_begin(_contact_name));
  const BoundaryElementMap::iterator bdend(env.boundary_elements_end(_contact_name));

  // min and max
  Point pmin(0);
  Point pmax(0);

  BoundaryElementMap::iterator bdit(bdfirst);
  for ( ; bdit != bdend; ++bdit)
  {
    const Elem* elem = *bdit;

  }

  SimulationEnvironment::ConstElemIterator it(env.elements_begin());
  SimulationEnvironment::ConstElemIterator end(env.elements_end());
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
  }

    // calculate min distance from contact
}



void
SchottkyTunneling::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();
  double ni = dd.get_intrinsic_density();
  double gn = 1; //dd.get_electron_gamma();
  double gp = 1; //dd.get_hole_gamma();

  recomb_e = recomb_h = (n * p - ni * ni * gn * gp);
}



void
SchottkyTunneling::get_net_recombination_rate_derivatives(
    std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  double n  = dd.get_electron_density();
  double p  = dd.get_hole_density();

  recomb_e[0] = recomb_h[0] =  p; // dR/dn
  recomb_e[1] = recomb_h[1] =  n; // dR/dp
}


/*
void
SchottkyTunneling::do_reinit(void)
{
  SimulationInterface* sim = SimulationInterface::get_simulation(get_simulator_id());
  unsigned int seqnum = sim->get_solve_sequence_number();

    {
      map<ID, vector<double> > data;
      data[_rec_id];

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
        AutoPtr<FEBase> fe(sim->build_finite_element(dim, FEType()));
        AutoPtr<QBase> qrule(QBase::build(libMeshEnums::QGAUSS, dim, libMeshEnums::FIFTH));
        fe->attach_quadrature_rule(qrule.get());

        const vector<Real>& JxW = fe->get_JxW();

        double tot_rec = 0.0;

        SimulationEnvironment::ConstElemIterator it(env.elements_begin());
        SimulationEnvironment::ConstElemIterator end(env.elements_end());
        for ( ; it != end; ++it)
        {
          const Elem* elem = *it;
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
*/
