/*  
 * This file is part of the tiberCAD module dd_generic.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file SchottkyTunneling.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "SchottkyTunneling.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/io/Database.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/io/Messages.h"
#include "tibercad/geom/MeshUtils.h"

#include "mesh_base.h"
#include "fe.h"
#include "quadrature_gauss.h"


#include <limits>

#include "tibercad/module/TiberModule.h"

//TIBER_MODULE(SchottkyTunneling, recombination, schottky_tunneling)

using namespace std;




SchottkyTunneling::SchottkyTunneling(const ModelOptions& options)
  : RecombinationModelInterface(options),
    _max_tunnel_length(20),
    _contact_voltage(0.0),
    _band('c'),
    _barrier(0.6),
    _mass(1.0)
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
  RecombinationModelInterface::do_init();

  _max_tunnel_length = get_option("maximum_tunnel_length", _max_tunnel_length);
  _contact_name = get_option("contact", _contact_name);
  get_parameter("effective_mass", _mass);

  if (_contact_name.empty())
    throw ModelErrorException("Need a contact name for Schottky "
        "contact tunneling model.");

  if (get_carrier_names().size() != 1)
    throw InitFailedException("Schottky tunneling model needs exactly "
        "one carrier specified");

  // calculate distance from contact and put it into
  // _elem_map if it is below the maximum tunnel length

  SimulationInterface* sim = SimulationInterface::get_simulation(get_simulator_id());
  const SimulationEnvironment& env = sim->get_environment();

  // the tunneling length in mesh units
  double tun_len = _max_tunnel_length * 1e-9 / sim->get_mesh_units();

  unsigned int dim = sim->get_mesh().mesh_dimension();

  std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, libMesh::FEType(libMesh::CONSTANT, libMesh::MONOMIAL)));
  libMesh::QGauss qrule(dim - 1, libMesh::CONSTANT);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Point>& normal = fe->get_normals();

  // first, build a bounding box for the contact to be faster afterwards
  SimulationEnvironment::BoundarySideIterator bdfirst(
      env.boundary_sides_begin(_contact_name));
  SimulationEnvironment::BoundarySideIterator bdend(
      env.boundary_sides_end(_contact_name));

  // min and max
  double max_double = numeric_limits<double>::max();
  Point pmin(max_double, max_double, max_double);
  Point pmax(-max_double, -max_double, -max_double);

  //Utils::Timer tt;

  // here we also extract the voltage option from the contact
  SimulationEnvironment::BoundarySideIterator bdit(bdfirst);
  for ( ; bdit != bdend; ++bdit)
  {
    const Elem* elem = (*bdit).elem();
    unsigned int side = (*bdit).side();

    if (bdit == bdfirst)
    {
      const PhysicalModel* mod =
          sim->get_interface_model<PhysicalModel>(elem, side);

      if (!has_option("voltage"))
        get_options().set_option("voltage",
            mod->get_options().get_option("voltage", string()));

      if (!has_option("barrier"))
      {
        if (mod->get_options().find_option("barrier_height"))
          get_options().set_option("barrier",
              mod->get_options().get_option("barrier_height", string()));
        if (mod->get_options().find_option("barrier"))
          get_options().set_option("barrier",
              mod->get_options().get_option("barrier", _barrier));
      }
      
      if (!has_option("band"))
      {
        string band = mod->get_options().get_option("band", "c");
        _band = band[0];
      }
      else
      {
        string band = get_option("barrier", "c");
        _band = band[0];
      }
    }

    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
      if (elem->is_node_on_side(i, side))
      {
        const Point& p = elem->point(i);
        pmin(0) = (p(0) - tun_len < pmin(0)) ? p(0) - tun_len : pmin(0);
        pmin(1) = (p(1) - tun_len < pmin(1)) ? p(1) - tun_len : pmin(1);
        pmin(2) = (p(2) - tun_len < pmin(2)) ? p(2) - tun_len : pmin(2);
        pmax(0) = (p(0) + tun_len > pmax(0)) ? p(0) + tun_len : pmax(0);
        pmax(1) = (p(1) + tun_len > pmax(1)) ? p(1) + tun_len : pmax(1);
        pmax(2) = (p(2) + tun_len > pmax(2)) ? p(2) + tun_len : pmax(2);
      }
  }

  //cerr << "Bounding box took " << tt.elapsed_string() << endl;
  //tt.reset();

  //cerr << "Bounding box: " << pmin << pmax << endl;

  SimulationEnvironment::ConstElemIterator it(env.elements_begin());
  SimulationEnvironment::ConstElemIterator end(env.elements_end());
  for ( ; it != end; ++it)
  {
    const Elem* elem = *it;
    const Point& centr = elem->vertex_average();

    // check if it is inside the tunneling bounding box
    if ((centr(0) < pmin(0)) || (centr(1) < pmin(1)) || (centr(2) < pmin(2)) ||
        (centr(0) > pmax(0)) || (centr(1) > pmax(1)) || (centr(2) > pmax(2)))
      continue;

    double mindist = numeric_limits<double>::max();
    Point p_mindist(mindist, mindist, mindist);

    // calculate min distance from contact
    for (bdit = bdfirst; bdit != bdend; ++bdit)
    {
      const Elem* elem = (*bdit).elem();
      unsigned int side = (*bdit).side();

      // TODO for now we use the distance between centroids, since we do not
      //      need the exact distance for the most primitive model
      Point normal(MeshUtils::get_outer_normal(elem, side));

      // side center as mean value of coordinates
      Point side_centr;
      int n_nodes = 0;
      for (int i = 0; i < elem->n_nodes(); ++i)
        if (elem->is_node_on_side(i, side))
        {
          side_centr += elem->point(i);
          ++n_nodes;
        }

      side_centr /= n_nodes;

      Point dist_vec(centr - side_centr);
      if ((dist_vec * normal) < 0)
      {
        double dist = dist_vec.norm();
        if (dist < mindist)
        {
          p_mindist = dist_vec;
          mindist = dist;
        }
      }
    }

    // NOTE: we insert the top parent, assuming that this list
    // is created before refinement
    _elem_map.insert(make_pair(elem->top_parent(), p_mindist));

  }

  //cerr << "Search took " << tt.elapsed_string() << endl;

  // now, extract the contact voltage we copied from the contact model
  get_parameter("voltage", _contact_voltage);

  // the same for the barrier height
  if (!has_option("barrier"))
    throw ModelErrorException("Need barrier height for Schottky "
        "contact tunneling model (contact: " + _contact_name + ")");
  _barrier = get_option("barrier", _barrier);

  // let the model know that it is a contact tunneling model
  set_tunneling_contact(env.get_boundary(_contact_name));

}



void
SchottkyTunneling::calculate_rate_and_derivatives(std::vector<double>& R,
    std::vector<std::vector<double>>& dPotentials)
{
  const ID id = this->get_carrier_ids()[0];

  const DriftDiffusionProperties& dd = get_driftdiffusionproperties();
  double Ec  = dd.get_carrier_properties(id)->get_band_edge();
  double Ef = -dd.get_q_fermi_potential(id);
  double kT = dd.get_lattice_temperature();
  //kT = dd.get_point_data().carrier_vt[0];

  double n  = dd.get_q_density(id);
  double dn  = dd.get_q_density_derivative(id);
  double qn = dd.get_carrier_properties(id)->get_charge();


  HashMap<const libMesh::Elem*, libMesh::Point>::Type::iterator it(
      _elem_map.find(dd.get_element()->top_parent()));

  // if the current element is not in our list, we can return immediately
  if (it == _elem_map.end()) return;

  double Gtun = 0.0;
  double band_edge = Ec - dd.get_electric_potential();;
  double pot_diff = 0.0;
  if (qn <= 0)
  {
    pot_diff = -_contact_voltage + _barrier - band_edge;
  }
  else
  {
    pot_diff =  _barrier + band_edge + _contact_voltage;
  }

  double E = dd.get_electric_field() * it->second / it->second.norm();
  E = 100 * E;

  if (pot_diff > 0.0)
  {
    double hcube = Constants::h * Constants::h * Constants::h;
    double A = 4 * M_PI * Constants::electron_mass * _mass *
        Constants::e * Constants::e * kT / hcube;

    double exp1, exp2;
    if (qn <= 0)
    {
      exp1 = exp(-(band_edge - Ef) / kT);
      exp2 = exp(-(band_edge + _contact_voltage) / kT);
    }
    else
    {
      exp1 = exp((band_edge - Ef) / kT);
      exp2 = exp((band_edge + _contact_voltage) / kT);

      E *= -1;
    }

    if (E >= 0) E = -1e-6;

    // E is negative !!
    double tmp = sqrt(2 * Constants::e * _mass * Constants::electron_mass *
        pot_diff * pot_diff * pot_diff) / E;
    double gamma = exp(4 / 3 / Constants::hbar * tmp);
    Gtun = -A * gamma * log((1 + exp1) / (1 + exp2)) * E / 1e6;

  }

  R[id] = Gtun;

}


