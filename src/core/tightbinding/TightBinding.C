//modules includes
#include "TightBinding.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModel.h"
#include "ElectricalContact.h"
#include "SimulationOptions.h"
#include "SimulationEnvironment.h"
#include "AtomisticStructure.h"


//libmesh includes
#include "mesh.h"

//-----------------------------------------------------------------------

TightBinding::TightBinding(){
  _atomistic_structure = NULL;
}


TightBinding::~TightBinding(){}


void 
TightBinding::do_init(){

  std::cerr << "Tight Binding Simulation Inizialization..." << std::endl;

  std::string fake_option;

  fake_option = get_options().get_option("fake_option","yo");
  std::cout << "fake_option is " << fake_option << std::endl;

}


void 
TightBinding::do_solve(){}


void 
TightBinding::parse_options(){}


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
 
  const std::string& modelname = options.get_option("type", "Heat_reservoir");
  
    ElectricalContact* model = ElectricalContact::create(modelname, options);
 
    if (model == NULL)  
      throw ModelErrorException("TightBinding: No such boundary model: " + modelname);
 
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
