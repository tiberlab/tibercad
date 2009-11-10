// $Id$

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

  destroy(_epsilon_model);

}

//=======================================================//
PhysicalModelInterface* MaxwellPhysicalModel::create_new (void) const
{
  return new MaxwellPhysicalModel();
}

//=======================================================//

void MaxwellPhysicalModel::do_init()
{
  if (_epsilon_model == NULL)
  {
    const ModelOptions& opt =  get_options ();

    _epsilon_model = OptDielectricConstant::create(get_material()->get_structure(), opt);

    _epsilon_model->set_material(get_material());

    _epsilon_model->init();
  }

}

//=======================================================//

void MaxwellPhysicalModel::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
  const MaxwellPhysicalModel* matA = dynamic_cast< const  MaxwellPhysicalModel*> (comp_A);
  const MaxwellPhysicalModel* matB = dynamic_cast< const  MaxwellPhysicalModel*> (comp_B);

  destroy(_epsilon_model);
  _epsilon_model = create_submodel_alloy(matA->_epsilon_model, matB->_epsilon_model, xa);
}
