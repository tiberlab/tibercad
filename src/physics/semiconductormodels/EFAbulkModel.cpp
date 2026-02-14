/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file EFAbulkModel.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/physics/semiconductormodels/EFAbulkModel.h"
#include "tibercad/physics/Material.h"


EFAbulkModel::EFAbulkModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}


//===================================//

void EFAbulkModel::do_print_info()
{
  //std::cout<<"(EFAbulkModel) "<< get_type() << "  "
  //         <<_bulkHamiltonian->get_type()<<" INFO:"<<std::endl;

  _bulkHamiltonian->print_info();
}


void EFAbulkModel::prepare_submodels(void)
{
  if (_bulkHamiltonian == nullptr)
  {
    const ModelOptions& opt =  get_options ();

    _bulkHamiltonian = EFAbulkHamiltonian::create(get_material(), opt);
    add_submodel("bulk_hamiltonian", _bulkHamiltonian);
  }
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



}

//====================================//

void EFAbulkModel::do_init_alloy (const PhysicalModel *comp_A, const PhysicalModel *comp_B, double xa)
{
  const EFAbulkModel* matA = dynamic_cast<const EFAbulkModel* > (comp_A);

  const EFAbulkModel* matB = dynamic_cast<const EFAbulkModel* > (comp_B);

  //destroy(_bulkHamiltonian);
  //_bulkHamiltonian = static_cast<EFAbulkHamiltonian*>(matA->_bulkHamiltonian->copy());
  //assert(_bulkHamiltonian != nullptr);
  //_bulkHamiltonian->set_owner(get_owner());
  //_bulkHamiltonian->init_alloy(matA->_bulkHamiltonian, matB->_bulkHamiltonian,xa);
}
