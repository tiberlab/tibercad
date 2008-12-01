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
TightBinding::build_rho(const double x, const double y, const double z)
{
  const double deltar_max = 7.5; //Maximum cutoff distance in Amstrong
  double deltar, uhatom;
  double rho = 0.0;
  double x1, y1, z1;


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



