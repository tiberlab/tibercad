/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file TightBinding.C
 * \brief tiberCAD API implementation.
 */



//modules includes
#include "tibercad/physics/tightbinding/TightBinding.h"
#include "tibercad/physics/BoundaryProperties.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/physics/tightbinding/TightBindingModelInterface.h"
#include "tibercad/base/SimulationOptions.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/atomistic/AtomisticStructure.h"
#include "tibercad/atomistic/BondMap.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/physics/PotentialInterface.h"
#include "tibercad/geom/QuantumContact.h"


//libmesh includes
#include "libmesh/mesh.h"

#include <map>

using namespace std;

//-----------------------------------------------------------------------

TightBinding::TightBinding(const ModelOptions& options)
 : EigenvalueProblem(options),
   _mulliken_netcharges()
{
  _pot_min = 0.0;
}


TightBinding::~TightBinding()
{
  _pot_shift.clear();
  _el_chem_pot.clear();
  _hl_chem_pot.clear();
  _mulliken_netcharges.clear();
  _u_hub.clear();
}


void
TightBinding::do_init()
{
  init_kspace(ModelOptions());
}


void
TightBinding::do_solve(){}


void
TightBinding::parse_options(){}


void
TightBinding::obtain_hubbard_parameters(){}





double
TightBinding::build_rho(const Elem* elem, const Point& r)
{
  AtomisticStructure *as = get_atomistic_structure();
  double scale = as->get_scale();
  const double deltar_max = 7.5 / scale; //Maximum cutoff distance in Amstrong
  double deltar, uhatom;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;

  const std::vector<Atom>& structure = as->get_structure_atoms();
  const std::vector<unsigned int>& atoms = get_elem_atoms(elem->id());

  x = r(0); y = r(1); z = r(2);

  if (_mulliken_netcharges.size() == 0)
    {
      std::cerr << "ERROR IN TIGHTBINDING: trying to build charge density "
      "but no mulliken charges are available" << std::endl;
      exit(1);
    }


  for (unsigned int id = 0; id  < atoms.size(); id++)
    {
      //Getting Hubbard parameter
      //Up to now densities are mapped on orbital S
      unsigned int iatm = atoms[id];
      Specie sp = structure[iatm].get_specie();
      uhatom = _u_hub[sp][S];

      //Convert atom position to mesh units
      x1 = structure[iatm].get_position(0) / scale;
      y1 = structure[iatm].get_position(1) / scale;
      z1 = structure[iatm].get_position(2) / scale;

      //delta_r is already in mesh units in this way
      deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));

      //Also hubbard parameters (and tau) must be scaled in mesh units
      // (uhatom is in (atomic units)^(-1))
      double tau = ( ( uhatom * ( 16.0 / 5.0 ) ) / (Constants::bohr_radius ) ) * get_mesh_units();

      if (deltar > deltar_max) continue;
      else
        {
	  rho = rho + (Constants::e * _mulliken_netcharges[iatm] * tau * tau * tau * exp(-1.0 * tau * deltar));

        }
    }

  rho = -rho / (8.0 * 3.141592653589793);

  //scale rho to C/cm^3
  double mesh_units = 100.0 * get_mesh_units();
  rho =  rho / (mesh_units * mesh_units * mesh_units);

  return rho;

}

void
TightBinding::get_solution_secure(const Elem* elem,
    const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
{

  std::vector<Point> points(elem->n_nodes());

  for (unsigned n = 0 ; n< elem->n_nodes(); ++n)
    {
      points[n] = elem->point(n);
    }

  get_solution_secure(elem,points,ids,values);

}


void
TightBinding::get_solution_secure(const Elem* elem, const std::vector<Point>& p,
    const std::set<ID>& ids, std::vector<std::map<ID, double> >& values)
{
  unsigned int np = p.size();
  values.resize(np);

  if (ids.count(CHARGE))
  {
    for (unsigned int n = 0; n < np; n++)
    {
      values[n][CHARGE] = build_rho(elem, p[n]);
    }
  }
  
}


void
TightBinding::project_potential(const std::string model_name, const std::string mode)
{
  Point p;

  // Use PotentialInterface to get the right simulation
  // The default is to assume them to com from the driftdiffusion module
  // (this is for backwards compatibility, but it would be better to always
  // specify the sources explicitly in the input)
  string el_pot  = get_option("electrostatic_potential", model_name);
  string el_chem = get_option("el_electrochemical_potential", model_name);
  string hl_chem = get_option("hl_electrochemical_potential", model_name);

  PotentialInterface pot_model(el_pot, "ElPotential");
  PotentialInterface elchem_model(el_chem, "eQFermi");
  PotentialInterface hlchem_model(hl_chem, "hQFermi");

  if( !pot_model.get_simulation()->is_solved() )
    throw InitFailedException("Potential model has not been solved");

  if (mode == "point")
  {//In point mode potential on atom is just kept as value on atom position
    //vector returned is sized number of atoms

    // atoms might be in quantum contacts, so we first obtain all IDs of
    // quantum contacts
    map<ID, QuantumContact*> quantum_contacts;
    Device& dev = get_environment().get_device();
    Device::quantum_contact_iterator qit(dev.quantum_contacts_begin());
    for ( ; qit != dev.quantum_contacts_end(); ++qit)
    {
      QuantumContact* qct = qit->second;
      quantum_contacts[qct->get_id()] = qct;
    }


    _pot_shift.clear();
    _pot_shift.resize(get_atomistic_structure()->get_N_atoms(), 0.0);
    _el_chem_pot.clear();
    _el_chem_pot.resize(get_atomistic_structure()->get_N_atoms(), 0.0);
    _hl_chem_pot.clear();
    _hl_chem_pot.resize(get_atomistic_structure()->get_N_atoms(), 0.0);

    // maybe we should take the mesh from the PotentialInterface?
    unsigned int dim = get_mesh().mesh_dimension();

    for (unsigned int i = 0; i < _pot_shift.size(); i++)
    {
      const Elem* elem = get_atomistic_structure()->get_structure_atom(i).get_elem();
      if (elem != NULL)
      {

        p(0) = get_atomistic_structure()->get_structure_atoms()[i].get_position(0)
                    / get_atomistic_structure()->get_scale();
        p(1) = get_atomistic_structure()->get_structure_atoms()[i].get_position(1)
		    / get_atomistic_structure()->get_scale();
        p(2) = get_atomistic_structure()->get_structure_atoms()[i].get_position(2)
		    / get_atomistic_structure()->get_scale();

        if (dim == 1) {p(1) = 0.0; p(2) = 0.0;}
        if (dim == 2) {p(2) = 0.0;}


        unsigned int subdomain = elem->subdomain_id();
        if (quantum_contacts.count(subdomain))
        {
          // this atom is in a quantum contact!
          QuantumContact* qct = quantum_contacts[subdomain];
          pair<const Elem*, Point > projected(qct->project_on_boundary(elem, p));
          elem = projected.first;
          p = projected.second;
        }

        // pot_shift is without "-" because the minus-sign is explicitly set in the
        // TB-codes.
        _pot_shift[i] = pot_model.get_potential(elem, p);
        _el_chem_pot[i] = elchem_model.get_potential(elem, p);
        _hl_chem_pot[i] = hlchem_model.get_potential(elem, p);

      }

    }
    //If atom has no element assigned, assigned the potential of the nearest neighbour
    //with non NULL element assigned

    for (unsigned int i = 0; i < _pot_shift.size(); i++)
    {
      //TODO: it works only for no-preserve and hydrogenation
      //we need to extend exploring neighbours until we don't reach one with non NULL
      //(it's not assured it will be the first one)
      if (get_atomistic_structure()->get_structure_atoms()[i].get_elem() == NULL)
      {
        int neighbour = get_atomistic_structure()->get_bond_map()[i][0];
        _pot_shift[i] = _pot_shift[neighbour];
        _el_chem_pot[i] = _el_chem_pot[neighbour];
        _hl_chem_pot[i] = _hl_chem_pot[neighbour];
      }
    }

  }

  //Process potential values to shift the smallest value to 0
  _pot_min = 0;

}




