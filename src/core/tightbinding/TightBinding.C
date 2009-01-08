//modules includes
#include "TightBinding.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModel.h"
#include "SimulationOptions.h"
#include "SimulationEnvironment.h"
#include "AtomisticStructure.h"


//libmesh includes
#include "mesh.h"

//-----------------------------------------------------------------------

TightBinding::TightBinding()
:_atomistic_structure(NULL),
_mulliken_netcharges()
{
}


TightBinding::~TightBinding()
{
}


void
TightBinding::do_init(){

  std::cerr << "Tight Binding Simulation Initialisation..." << std::endl;

}


void
TightBinding::do_solve(){}


void
TightBinding::parse_options(){}

void
TightBinding::obtain_hubbard_parameters(){}


PhysicalModel*
TightBinding::create_physical_model (const ModelOptions &options,
    const Material* mat) const throw (ModelErrorException)
    {

      TightBindingModel* model = dynamic_cast<TightBindingModel*> ( PhysicalModelInterface::create("tightbinding",options) );

      if (model == NULL)
        throw ModelErrorException("TightBinding: Tight Binding physical model is not created" );

      return model;

    }



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

  AtomisticStructure* atomistic_structure = NULL;

  if (get_options().find_option("atomistic_structure") )
    {
      std::string name;
      name = get_options().get_option("atomistic_structure", "none");
      if (name.compare("none") != 0){
        _atomistic_structure = get_environment().get_device().get_atomistic_structure(name);
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
  const double deltar_max = 7.5; //Maximum cutoff distance in Amstrong
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

      std::cout << "rho before loop is " << rho << std::endl;
      //Getting Hubbard parameter
      //Up to now densities are mapped on orbital S
      uhatom = _u_hub[_atomistic_structure->get_structure_atoms()[iatm].get_specie()][S];

      x1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(1);
      y1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(2);
      z1 = _atomistic_structure->get_structure_atoms()[iatm].get_position(3);

      deltar = sqrt( (x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));

      std::cout << "deltar is " << deltar <<std::endl;

      if (deltar > deltar_max) continue;
      else
        {
          std::cout << "uhatom is " << uhatom <<std::endl;
          rho = rho + 16.384 * _mulliken_netcharges[iatm] * uhatom * uhatom * uhatom * exp(-3.20*uhatom*deltar);
        }
      std::cout << "rho after loop is " << rho << std::endl;
    }
  rho = -rho / 4.0 / 3.141592653589793;

  std::cout << "rho finally is " << rho << std::endl;

  return rho;

}


ID
TightBinding::convert_variable_name_to_id(const std:: string& variable_name) const
{

  ID id = INVALID_ID;


  if (variable_name == "charge_density" )
    id  = CHARGE;


  return id;
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

  if (ids.count(CHARGE))
    {
      for (unsigned int n = 0; n < p.size(); n++)
        values[n][CHARGE] = build_rho(p[n]);
    }

}




void
TightBinding::build_elemental_results(const std::set<std::string>& variables,
    std::vector<double>& results, std::vector<std::string>& legend)
{


  // we only do something if we are on processor 0
  // TODO parallelize
  if (libMesh::processor_id() != 0)
    return;

  unsigned int n_vars = 0;
  const unsigned int nn  = _mesh->n_active_elem();
  int ch = -1;

  if (variables.count("QuantumCharge"))
    {
      legend.resize( 1 );
      legend[0] = "qDensity";
      ch = n_vars;
      n_vars++;
    }


  legend.resize(n_vars);

  results.resize(nn * n_vars,0.0);

  MeshBase::const_element_iterator it =  _mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end =     _mesh->active_local_elements_end();


  unsigned int elem_number = 0;
  for ( ; it != end; ++it)
    {

      const Elem* elem = *it;

      unsigned int id = n_vars * elem_number;

      std::vector<std::map<ID, double> > values;
      std::set<ID> ids;
      ids.insert(id);

      get_solution_secure(elem, ids, values);

      double charge = values[0][CHARGE];


      if (ch != -1)
        {
          results[id + ch] = charge;
        }


      elem_number++;
    } //over element

  results.resize(elem_number * n_vars);

}
