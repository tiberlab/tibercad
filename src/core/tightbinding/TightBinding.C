

//modules includes
#include "TightBinding.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "TightBindingModel.h"
#include "ElectricalContact.h"
#include "SimulationOptions.h"
#include "DftbpWrapper.h"


//libmesh includes
#include "mesh.h"

using namespace std;

//-----------------------------------------------------------------------

TightBinding::TightBinding(){}

TightBinding::~TightBinding(){}


void TightBinding::do_init(){
  cout << "TIGHT BINDING INITIALIZATION, I'll create a DFTBp instance " << endl;
  DftbpWrapper*  inst1;
  //DFTBp_instance = new(DftbpWrapper);
  inst1 = DftbpWrapper::create();

 int nAtom = 2;
  double eTemp = 1e-8;
  int iPeriodic = 1;
  int nType = 1;
  char *speciesNames;
  int *species;
  
  speciesNames = (char *) malloc(nType * DFTBP_MC * sizeof(char)); 
  strcpy(&speciesNames[0], "Si");
  for (int ii = strlen(&speciesNames[0]); ii < DFTBP_MC; ++ii) {
    speciesNames[0+ii] = DFTBP_PADCHAR;
  }

  species = new int[nAtom];
  species[0] = 1;
  species[1] = 1; 

  cout << "Calling fill_param... ";
  inst1->fill_param(nAtom, nType, eTemp, iPeriodic, speciesNames, species);
  cout << "done." << endl;
  delete species;
  free(speciesNames);

  char *skNames;
  skNames = (char *) malloc(1 * DFTBP_LC * sizeof(char));
  strcpy(skNames, "~/Si-Si.skf");
  for (int ii = strlen(skNames); ii < DFTBP_LC; ++ii) {
    skNames[ii] = DFTBP_PADCHAR;
  }

  int *mAngs;
  mAngs = new int[nType];
  mAngs[0] = 1;
  int orbResolved = 0;
  int skInterp = 2;
  cout << "Calling addskdata... ";
  inst1->addskdata(skNames, mAngs, orbResolved, skInterp, nType);
  cout << "done." << endl;
  delete mAngs;
  free(skNames);

  double *latVecs; 
  latVecs = new double[9]; 
  latVecs[0] = 5.12785877;
  latVecs[1] = 5.12785877;
  latVecs[2] = 0.00000000;
  latVecs[3] = 0.00000000;
  latVecs[4] = 5.12785877;
  latVecs[5] = 5.12785877;
  latVecs[6] = 5.12785877;
  latVecs[7] = 0.00000000;
  latVecs[8] = 5.12785877;

  cout << "Calling addlattice... ";
  inst1->addlattice(latVecs);
  cout << "done." << endl;
  delete latVecs;

  int nKPoint = 4;
  double *kPoints;
  kPoints = new double[3*nKPoint];
  kPoints[0] = 0.25;
  kPoints[1] = 0.25;
  kPoints[2] = 0.25;
  kPoints[3] = -0.25;
  kPoints[4] = 0.25;
  kPoints[5] = 0.25;
  kPoints[6] = 0.25;
  kPoints[7] = -0.25;
  kPoints[8] = 0.25;
  kPoints[9] = -0.25;
  kPoints[10] = -0.25;
  kPoints[11] = 0.25;

  double *kWeights;
  kWeights = new double[nKPoint];
  kWeights[0] = 1.0;
  kWeights[1] = 1.0;
  kWeights[2] = 1.0;
  kWeights[3] = 1.0;

  cout << "Calling addkpoints... ";
  inst1->addkpoints(nKPoint, kPoints, kWeights);
  cout << "done." << endl;
  delete kPoints;
  delete kWeights;

  cout << "Feeding DFTB with the input data... ";
  inst1->initdftb();
  cout << "done." << endl;

  double *coords;
  coords = new double[3*nAtom];
  coords[0] = 0.0;
  coords[1] = 0.0;
  coords[2] = 0.0;
  coords[3] = 2.563929;
  coords[4] = 2.563929;
  coords[5] = 2.563929;

  cout << "Calling up_coords... ";
  inst1->up_coords(nAtom, coords);
  cout << "done." << endl;
  delete coords;

  double energy = 0;
  cout << "Calling get_energy... ";
  inst1->get_energy(energy);
  cout <<  "done." << endl;

  cout << endl << "DFTB ran correctly: energy computed is " << energy << endl
       << endl;

  cout << "Now I destroy it" << endl;
  delete inst1;
  cout << "Ok" << endl;
  cout << endl << endl << endl;
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
