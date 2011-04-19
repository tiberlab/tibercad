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


//libmesh includes
#include "mesh.h"

//-----------------------------------------------------------------------

TightBinding::TightBinding(const ModelOptions& options)
 : EigenvalueProblem(options),
   _atomistic_structure(NULL),
   _mulliken_netcharges(),
   _mesh(NULL)
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
  //Get mesh reference
  _mesh = & ( get_environment().get_device().get_mesh());

  // Getting reference to atomistic structure for calculation
  get_atomistic_structure();
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


void
TightBinding::get_atomistic_structure(void){

  //Make a local copy of atomistic structure, in order to perform operations as
  //strain dependent atom displacement
  AtomisticStructure* atomistic_structure = NULL;

  if (get_options().find_option("atomistic_structure") )
    {
      std::string name;
      name = get_options().get_option("atomistic_structure", "none");
      if (name.compare("none") != 0){
       // _atomistic_structure = new AtomisticStructure( *(get_environment().get_device().get_atomistic_structure(name)) );
        _atomistic_structure = (get_environment().get_device().get_atomistic_structure(name));
      }
    }
  else
    {
      std::cerr << "ERROR in Tight Binding Simulation: an atomistic structure  must be specified "
      << get_name() << std::endl;
      exit(0);
    }
}


double
TightBinding::build_rho(const Point& r)
{
  const double deltar_max = 7.5 / _atomistic_structure->get_scale(); //Maximum cutoff distance in Amstrong
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

  for (unsigned int iatm = 0; iatm  < _atomistic_structure->get_N_atoms(); iatm++)
    {
      //Getting Hubbard parameter
      //Up to now densities are mapped on orbital S
      Specie sp = _atomistic_structure->get_structure_atoms()[iatm].get_specie();
      uhatom = _u_hub[sp][S];

      //Convert atom position to mesh units
      x1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(1) / _atomistic_structure->get_scale();
      y1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(2) / _atomistic_structure->get_scale();
      z1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(3) / _atomistic_structure->get_scale();

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

  std::cout << "setting simulation " << std::endl;
  model.set_simulation(model_name);

  if( !model.get_simulation()->is_solved() )
    throw InitFailedException("Potential model has not been solved");

  if (mode == "point")
    {//In point mode potential on atom is just kept as value on atom position
      //vector returned is sized number of atoms

      _pot_shift.clear();
      _pot_shift.resize(_atomistic_structure->get_N_atoms(), 0.0);
      _el_chem_pot.clear();
      _el_chem_pot.resize(_atomistic_structure->get_N_atoms(), 0.0);
      _hl_chem_pot.clear();
      _hl_chem_pot.resize(_atomistic_structure->get_N_atoms(), 0.0);

      unsigned int dim = get_environment().get_device().get_mesh().mesh_dimension();

      for (unsigned int i = 0; i < _pot_shift.size(); i++)
	{
	  if (_atomistic_structure->get_structure_atoms()[i].get_elem() != NULL)
	    {
	      p(0) = _atomistic_structure->get_structure_atoms()[i].get_position()(1)
                / _atomistic_structure->get_scale();
	      p(1) = _atomistic_structure->get_structure_atoms()[i].get_position()(2)
		/ _atomistic_structure->get_scale();
	      p(2) = _atomistic_structure->get_structure_atoms()[i].get_position()(3)
		/ _atomistic_structure->get_scale();

	      if (dim == 1) {p(1) = 0.0; p(2) = 0.0;}
	      if (dim == 2) {p(2) = 0.0;}


	      // pot_shift is without "-" because the minus-sign is explicitly set in the
	      // TB-codes. 
	      _pot_shift[i] = model.get_potential(_atomistic_structure->
						  get_structure_atoms()[i].get_elem(), p);
	      _el_chem_pot[i] = model.get_el_chem_potential(_atomistic_structure->
	                                          get_structure_atoms()[i].get_elem(), p);
	      _hl_chem_pot[i] = model.get_hl_chem_potential(_atomistic_structure->
	                                          get_structure_atoms()[i].get_elem(), p);
	      //std::cout << " shifting " << _pot_shift[i] << std::endl;
	    }
	  else
	    {
	      _pot_shift[i] = 0.0;
	      _el_chem_pot[i] = 0.0;
	      _hl_chem_pot[i] = 0.0;
	    }
	}
      //If atom has no element assigned, assigned the potential of the nearest neighbour
      //with non NULL element assigned

    for (unsigned int i = 0; i < _pot_shift.size(); i++)
      {
        //TODO: it works only for no-preserve and hydrogenation
        //we need to extend exploring neighbours until we don't reach one with non NULL
        //(it's not assured it will be the first one)
        if (_atomistic_structure->get_structure_atoms()[i].get_elem() == NULL)
          {
            int neighbour = _atomistic_structure->get_bond_map()[i][0];
            _pot_shift[i] = _pot_shift[neighbour];
            _el_chem_pot[i] = _el_chem_pot[neighbour];
            _hl_chem_pot[i] = _hl_chem_pot[neighbour];
          }
      }

    }

  //Process potential values to shift the smallest value to 0
  _pot_min = _pot_shift[0];
  double tmp;
  for (unsigned int i = 1; i < _pot_shift.size(); i++)
  {
  //TODO: change in a cycle without passivation atoms (and move to empirical tight binding module!),
  //comparison with "H" in any cycle is not a good idea
      if ( (_pot_shift[i] < _pot_min) && 
      (_atomistic_structure->get_structure_atoms()[i].get_specie() != "H")) _pot_min = _pot_shift[i];
  }
  //std::cout << "pot_min is " << _pot_min << std::endl;
  for (unsigned int i = 0; i < _pot_shift.size(); i++)
  {
      tmp = _pot_shift[i] - _pot_min;
      _pot_shift[i] = tmp;
      //Switch off potential projection for passivation atoms!!!!!! Potential shift is set at 0.0
      if (_atomistic_structure->get_structure_atoms()[i].get_specie() == "H")
      _pot_shift[i] = 0.0;
  }

  double* pot = new double[_pot_shift.size()];
      for (unsigned int i = 0; i < _pot_shift.size(); i++)
        pot[i] = _pot_shift[i];
      _atomistic_structure->print_structure("pot_on_atom.xyz", pot);

  std::cout << "project_potential done " << std::endl;
}




