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

using namespace std;

//-----------------------------------------------------------------------

TightBinding::TightBinding(){}

TightBinding::~TightBinding(){}


void TightBinding::do_init(){

  cerr << "Tight Binding Simulation Inizialization..." << endl;

  string fake_option;

  fake_option = get_options().get_option("fake_option","yo");
  cout << "fake_option is " << fake_option << endl;

}


void TightBinding::do_solve(){}


void TightBinding::parse_options(){}


PhysicalModel*   TightBinding::create_physical_model (const ModelOptions &options) const 
                    throw (ModelErrorException)
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
 
    const string& modelname = options.get_option("type", "Heat_reservoir");
  
    ElectricalContact* model = ElectricalContact::create(modelname, options);
 
    if (model == NULL)  
      throw ModelErrorException("TightBinding: No such boundary model: " + modelname);
 
   return model;
 
}
