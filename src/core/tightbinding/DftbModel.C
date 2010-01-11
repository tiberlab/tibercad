

//modules includes
#include "DftbModel.h"
#include "SimulationInterface.h"

using namespace std;

//-----------------------------------------------------------


DftbModel::DftbModel(const ModelOptions& options)
 : TightBindingModelInterface(options)
{
  cout << "DFTB Physical Model has been created" << endl;
}


DftbModel::~DftbModel(){}


void DftbModel::do_init(){}
