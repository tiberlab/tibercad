

//modules includes
#include "EtbModel.h"
#include "SimulationInterface.h"

using namespace std;

//-----------------------------------------------------------


ETBModel::ETBModel(const ModelOptions& options)
 : TightBindingModelInterface(options)
{
  cout << "Empirical Tight Binding Model has been created" << endl;
}


ETBModel::~ETBModel(){}


void ETBModel::do_init(){}
