#include "EFAbulkModel.h"
#include "Material.h"


EFAbulkModel::EFAbulkModel(const ModelOptions& options) :
  PhysicalModel(options)
{
  _bulkHamiltonian = NULL;
}

//===================================//

EFAbulkModel::~EFAbulkModel()
{
  destroy(_bulkHamiltonian);
}

//===================================//

void EFAbulkModel::do_print_info()
{
  //std::cout<<"(EFAbulkModel) "<< get_type() << "  "
  //         <<_bulkHamiltonian->get_type()<<" INFO:"<<std::endl;

  _bulkHamiltonian->print_info();
}

//===================================//
void EFAbulkModel::do_init()
{

  get_option("model","");
  get_option("particle","");
  get_option("kpVVtermSymmetric","");
  get_option("kpCVtermSymmetric","");
  get_option("temperature_scaling","");
  get_option("consider_temperature","");
  get_option("spurious","");

  if (_bulkHamiltonian == NULL)
  {
    const ModelOptions& opt =  get_options ();

    _bulkHamiltonian = EFAbulkHamiltonian::create(get_material(), opt);
    _bulkHamiltonian->set_owner(get_owner());
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
  _bulkHamiltonian->set_owner(get_owner());
  _bulkHamiltonian->init_alloy(matA->_bulkHamiltonian, matB->_bulkHamiltonian,xa);
}
