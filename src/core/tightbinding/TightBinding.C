// $Id$


//modules includes
#include "TightBinding.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModelInterface.h"
#include "SimulationOptions.h"
#include "SimulationEnvironment.h"
#include "AtomisticStructure.h"
#include "Constants.h"
#include "PotentialInterface.h"
#include "QuantumContact.h"


//libmesh includes
#include "mesh.h"

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
  init_kspace();
}


void
TightBinding::do_solve(){}


void
TightBinding::parse_options(){}


void
TightBinding::obtain_hubbard_parameters(){}



//Careful!!!! What means boundary condition in tight binding????????
BoundaryProperties* TightBinding::create_boundary_model (const ModelOptions &options) const
throw (ModelErrorException)

{
  BoundaryProperties* model = NULL;

  // these come from other models!
  /*
  const std::string& modelname = options.get_option("type", "Heat_reservoir");

    model = ElectricalContact::create(modelname, options);

    if (model == NULL)
      throw ModelErrorException("TightBinding: No such boundary model: " + modelname);
   */
  return model;

}




double
TightBinding::build_rho(const Point& r)
{
  const double deltar_max = 7.5 / get_atomistic_structure()->get_scale(); //Maximum cutoff distance in Amstrong
  double deltar, uhatom;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;

  x = r(0); y = r(1); z = r(2);

  if (_mulliken_netcharges.size() == 0)
    {
      std::cerr << "ERROR IN TIGHTBINDING: trying to build charge density "
      "but no mulliken charges are available" << std::endl;
      exit(1);
    }

  for (unsigned int iatm = 0; iatm  < get_atomistic_structure()->get_N_atoms(); iatm++)
    {
      //Getting Hubbard parameter
      //Up to now densities are mapped on orbital S
      Specie sp = get_atomistic_structure()->get_structure_atoms()[iatm].get_specie();
      uhatom = _u_hub[sp][S];

      //Convert atom position to mesh units
      x1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(0) / get_atomistic_structure()->get_scale();
      y1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(1) / get_atomistic_structure()->get_scale();
      z1 = get_atomistic_structure()->get_structure_atoms()[iatm].get_position(2) / get_atomistic_structure()->get_scale();

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


//ID
//TightBinding::convert_variable_name_to_id(const std:: string& variable_name) const
//{
//
//  ID id = INVALID_ID;
//
//
//  if (variable_name == "charge_density" )
//    id  = CHARGE;
//
//  return id;
//}


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
        values[n][CHARGE] = build_rho(p[n]);
        }
    }

}


//void
//TightBinding::build_elemental_results(const std::set<std::string>& variables,
//    std::vector<double>& results, std::vector<std::string>& legend)
//{
//
//
//  // we only do something if we are on processor 0
//  // TODO parallelize
//  if (libMesh::processor_id() != 0)
//    return;
//
//
//  // if there is no mesh we can return immediately
//  if (_mesh == NULL)
//    return;
//
//
//
//  unsigned int n_vars = 0;
//  const unsigned int nn  = _mesh->n_active_elem();
//  int ch = -1;
//
//  if (variables.count("QuantumCharge"))
//    {
//      legend.resize( 1 );
//      legend[0] = "qDensity";
//      ch = n_vars;
//      n_vars++;
//    }
//
//
//  legend.resize(n_vars);
//
//  results.resize(nn * n_vars,0.0);
//
//  MeshBase::const_element_iterator it =  _mesh->active_local_elements_begin();
//  const MeshBase::const_element_iterator end =     _mesh->active_local_elements_end();
//
//
//  unsigned int elem_number = 0;
//  for ( ; it != end; ++it)
//    {
//
//      const Elem* elem = *it;
//
//      unsigned int id = n_vars * elem_number;
//
//      std::vector<std::map<ID, double> > values;
//
//      std::set<ID> ids;
//      ids.insert(CHARGE);
//
//      get_solution_secure(elem, ids, values);
//
//      double charge = values[0][CHARGE];
//
//
//      if (ch != -1)
//        {
//          results[id + ch] = charge;
//        }
//
//
//      elem_number++;
//    } //over element
//
//  results.resize(elem_number * n_vars);
//
//}


void
TightBinding::project_potential(const std::string model_name, const std::string mode)
{
  Point p;

  //Use PotentialInterface to get the right simulation
  PotentialInterface model;

  model.set_simulation(model_name);

  if( !model.get_simulation()->is_solved() )
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
          pair<const Elem*, Point> projected(qct->project_on_boundary(elem, p));
          elem = projected.first;
          p = projected.second;
        }

        // pot_shift is without "-" because the minus-sign is explicitly set in the
        // TB-codes.
        _pot_shift[i] = model.get_potential(elem, p);
        _el_chem_pot[i] = model.get_el_chem_potential(elem, p);
        _hl_chem_pot[i] = model.get_hl_chem_potential(elem, p);
        //std::cout << " shifting " << _pot_shift[i] << std::endl;
      }
      //else
      //  {
      //    _pot_shift[i] = 0.0;
      //    _el_chem_pot[i] = 0.0;
      //    _hl_chem_pot[i] = 0.0;
      //  }
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

  /*
  _pot_min = _pot_shift[0];
  double tmp;
  for (unsigned int i = 1; i < _pot_shift.size(); i++)
  {
  //TODO: change in a cycle without passivation atoms (and move to empirical tight binding module!),
  //comparison with "H" in any cycle is not a good idea
      if ( (_pot_shift[i] < _pot_min) && 
      (get_atomistic_structure()->get_structure_atoms()[i].get_specie() != "H")) _pot_min = _pot_shift[i];
  }
  //std::cout << "pot_min is " << _pot_min << std::endl;
  for (unsigned int i = 0; i < _pot_shift.size(); i++)
  {
      tmp = _pot_shift[i] - _pot_min;
      _pot_shift[i] = tmp;
      //Switch off potential projection for passivation atoms!!!!!! Potential shift is set at 0.0
      if (get_atomistic_structure()->get_structure_atoms()[i].get_specie() == "H")
      _pot_shift[i] = 0.0;
  }

  double* pot = new double[_pot_shift.size()];
  for (unsigned int i = 0; i < _pot_shift.size(); i++)
    pot[i] = _pot_shift[i];
  get_atomistic_structure()->print_structure("pot_on_atom.xyz", pot);
  delete pot;
   */

}




void
TightBinding::build_map_elem_atoms(double projection_length)
{
  
  // Get total number of elements
  // (the map is oversized, but faster since the elem ID is used as key)
  _elem_to_atoms.resize(get_mesh().n_elem());

  double scale = get_atomistic_structure()->get_scale();

  // Maximum cutoff distance
  //const double tau = 1.0 / projection_length; // projection in Angstroms
  //const double deltar_max = (5.0*log(10.0) - log( 8.0*3.141593/(tau*tau*tau) )  ) / tau;   

  
  const double sigma = projection_length;
  const double sigma2 = 2.0*sigma*sigma;
  const double deltar2_max = sigma2*(5.0*log(10.0) - 1.5*log(2.0*3.141593*sigma)); 
  const double deltar_max = sqrt(deltar2_max);

  double deltar, deltar2, uhatom;
  double rho = 0.0;
  double x1, y1, z1;
  double x ,y, z;
  
  unsigned int N_wo_H = get_atomistic_structure()->get_N_without_H();

  // Estimate number of atoms in a sphere
  unsigned int Nat = round( sqrt(3.0)*3.1416/2.0 * pow(deltar_max/1.90, 3.0) );

  std::cout<<"(TB) proj length: "<<projection_length<<std::endl;
  std::cout<<"(TB) Rmax: "<<deltar_max<<std::endl;
  std::cout<<"(TB) Nat: "<<Nat<<std::endl;

  const std::vector<Atom>& structure = get_atomistic_structure()->get_structure_atoms();

  MeshBase::const_element_iterator it = get_mesh().elements_begin();
  const MeshBase::const_element_iterator end = get_mesh().elements_end();

  for ( ; it != end; ++it)
  { 
    const Elem* elem = *it;
    int id = elem->id();
    std::vector<unsigned int> temp;
    unsigned int count = 0;

    temp.reserve(Nat);
 
    Point p = elem->centroid();
    x=p(0)*scale;
    y=p(1)*scale;
    z=p(2)*scale;

    for (unsigned int iatm = 0; iatm  <  N_wo_H; iatm++)
    {
   
      x1 = structure[iatm].get_position(0);
      y1 = structure[iatm].get_position(1);
      z1 = structure[iatm].get_position(2);

      //deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));
      deltar2 = (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1);

      if (deltar2 > deltar2_max) continue;
      else
      {
        temp.push_back(iatm);
	count++;
      }
    }

    _elem_to_atoms[id].resize(count);

    for (unsigned int iatm = 0; iatm  <  count; iatm++)
    	_elem_to_atoms[id][iatm]=temp[iatm];

    temp.clear();
    
  }
}
