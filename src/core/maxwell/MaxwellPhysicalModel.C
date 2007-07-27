#include "MaxwellPhysicalModel.h"
#include "Material.h"

using namespace std;

//======================================================//
MaxwellPhysicalModel::MaxwellPhysicalModel()
{
  _epsilon_model = NULL;
}

//=======================================================//

MaxwellPhysicalModel::~MaxwellPhysicalModel()
{

  PhysicalModelInterface::destroy(_epsilon_model);

}

//=======================================================//
PhysicalModelInterface* MaxwellPhysicalModel::create_new (void) const
{
  return new MaxwellPhysicalModel();
}

//=======================================================//

void MaxwellPhysicalModel::do_init()
{
  PhysicalModelInterface::destroy(_epsilon_model);

  const ModelOptions& opt =  get_options ();

  _epsilon_model = OptDielectricConstant::create(get_material()->get_structure(), opt);

  _epsilon_model->init();

}

//=======================================================//

void MaxwellPhysicalModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
  const MaxwellPhysicalModel* matA = dynamic_cast< const  MaxwellPhysicalModel*> (comp_A);
  const MaxwellPhysicalModel* matB = dynamic_cast< const  MaxwellPhysicalModel*> (comp_B);
 
  
 
  _epsilon_model->build_alloy(matA->_epsilon_model, matB->_epsilon_model, xa);
}
