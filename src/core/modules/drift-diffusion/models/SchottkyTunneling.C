// $Id$

#include "SchottkyTunneling.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "mesh_base.h"
#include "fe.h"
#include "quadrature_gauss.h"

#include <limits>


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
  _max_tunnel_length = get_option("maximum_tunnel_length", _max_tunnel_length);
  _contact_name = get_option("contact", _contact_name);

  if (_contact_name.empty())
    throw ModelErrorException("Need a contact name for Schottky "
        "contact tunneling model.");

  // calculate distance from contact and put it into
  // _elem_map if it is below the maximum tunnel length

  SimulationInterface* sim = SimulationInterface::get_simulation(get_simulator_id());
  const SimulationEnvironment& env = sim->get_environment();

  // the tunneling length in mesh units
  double tun_len = _max_tunnel_length * 1e-9 / sim->get_mesh_units();

  unsigned int dim = sim->get_mesh().mesh_dimension();

  AutoPtr<FEBase> fe(FEBase::build(dim, FEType()));
  QGauss qrule(dim - 1, CONSTANT);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Point>& normal = fe->get_normals();

  // first, build a bounding box for the contact to be faster afterwards
  SimulationEnvironment::BoundarySideIterator bdfirst(
      env.boundary_sides_begin(_contact_name));
  SimulationEnvironment::BoundarySideIterator bdend(
      env.boundary_sides_end(_contact_name));

  // min and max
  Point pmin(numeric_limits<double>::max());
  Point pmax(-numeric_limits<double>::max());

  SimulationEnvironment::BoundarySideIterator bdit(bdfirst);
  for ( ; bdit != bdend; ++bdit)
  {
    const Elem* elem = (*bdit).elem();
    unsigned int side = (*bdit).side();

    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
      if (elem->is_node_on_side(i, side))
      {
        const Point& p = elem->point(i);
        //cerr << p << endl;
        pmin(0) = (p(0) - tun_len < pmin(0)) ? p(0) - tun_len : pmin(0);
        pmin(1) = (p(1) - tun_len < pmin(1)) ? p(1) - tun_len : pmin(1);
        pmin(2) = (p(2) - tun_len < pmin(2)) ? p(2) - tun_len : pmin(2);
        pmax(0) = (p(0) + tun_len > pmax(0)) ? p(0) + tun_len : pmax(0);
        pmax(1) = (p(1) + tun_len > pmax(1)) ? p(1) + tun_len : pmax(1);
        pmax(2) = (p(2) + tun_len > pmax(2)) ? p(2) + tun_len : pmax(2);
      }
  }

  //cerr << "Bounding box: " << pmin << pmax << endl;

  SimulationEnvironment::ConstElemIterator it(env.elements_begin());
  SimulationEnvironment::ConstElemIterator end(env.elements_end());
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    const Point& centr = elem->centroid();

    double mindist = numeric_limits<double>::max();
    Point p_mindist(mindist);

    // check if it is inside the tunneling bounding box
    if ((centr(0) < pmin(0)) || (centr(1) < pmin(1)) || (centr(1) < pmin(1)) ||
        (centr(0) > pmax(0)) || (centr(0) > pmax(0)) || (centr(0) > pmax(0)))
      continue;

    // calculate min distance from contact
    for (bdit = bdfirst; bdit != bdend; ++bdit)
    {
      const Elem* elem = (*bdit).elem();
      unsigned int side = (*bdit).side();

      fe->reinit(elem, side);
      // TODO for now we use the distance between centroids, since we do not
      //      need the exact distance for the most primitive model
      AutoPtr<Elem> side_el(elem->build_side(side));
      const Point& side_centr = side_el->centroid();

      Point dist(centr - side_centr);
      if ((dist * normal[0]) < 0)
        if (dist.size() < mindist)
          p_mindist = dist;
    }

    // NOTE: we insert the top parent, assuming that this list
    // is created before refinement
    _elem_map.insert(make_pair(elem->top_parent(), p_mindist));

  }

}



void
SchottkyTunneling::get_net_recombination_rates(double& recomb_e,
    double& recomb_h)
{
  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();

  recomb_e = recomb_h = 0.0;

  // if the current element is not in our list, we can return immediately
  if (!_elem_map.count(dd.get_element()->top_parent())) return;


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

  recomb_e[0] = recomb_h[0] =  0;
  recomb_e[1] = recomb_h[1] =  0;

  // if the current element is not in our list, we can return immediately
  if (!_elem_map.count(dd.get_element()->top_parent())) return;

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
