#include "EFAbulkModel.h"
#include "Material.h"


EFAbulkModel::EFAbulkModel()
{
  _bulkHamiltonian = NULL;
}

//===================================//

EFAbulkModel::~EFAbulkModel()
{
  destroy(_bulkHamiltonian);
}

//===================================//

void EFAbulkModel::do_init()
{

  if (_bulkHamiltonian == NULL)
  {
    const ModelOptions& opt =  get_options ();

    _bulkHamiltonian = EFAbulkHamiltonian::create(get_material() -> get_structure(), opt);

    _bulkHamiltonian->set_material(get_material());

    _bulkHamiltonian->init();
  }


}

//====================================//

void EFAbulkModel::do_init_alloy (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
  const EFAbulkModel* matA = dynamic_cast<const EFAbulkModel* > (comp_A);

  const EFAbulkModel* matB = dynamic_cast<const EFAbulkModel* > (comp_B);

  destroy(_bulkHamiltonian);
  _bulkHamiltonian = static_cast<EFAbulkHamiltonian*>(matA->_bulkHamiltonian->copy());
  assert(_bulkHamiltonian != NULL);
  _bulkHamiltonian->set_material(get_material());
  _bulkHamiltonian->init_alloy(matA->_bulkHamiltonian, matB->_bulkHamiltonian,xa);
}
