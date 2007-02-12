#include "EFAbulkModel.h"
#include "Material.h"


EFAbulkModel::EFAbulkModel()
{
  _bulkHamiltonian = NULL;
}

//===================================//

EFAbulkModel::~EFAbulkModel()
{
  PhysicalModelInterface::destroy(_bulkHamiltonian);
}

//===================================//

void EFAbulkModel::do_init()
{
  PhysicalModelInterface::destroy(_bulkHamiltonian);

  const ModelOptions& opt =  get_options ();

  _bulkHamiltonian = EFAbulkHamiltonian::create(get_material() -> get_structure(), opt);

  _bulkHamiltonian->set_material(get_material());

  _bulkHamiltonian->init();
  
  
}

//====================================//

void EFAbulkModel::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
  const EFAbulkModel* matA = dynamic_cast<const EFAbulkModel* > (comp_A);

  const EFAbulkModel* matB = dynamic_cast<const EFAbulkModel* > (comp_B);

  _bulkHamiltonian->build_alloy(matA->_bulkHamiltonian, matB->_bulkHamiltonian,xa);
}
