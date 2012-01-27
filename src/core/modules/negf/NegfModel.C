/*
 * Negfmodel.C
 *
 *  Created on: Jan 18, 2012
 *      Author: fpalomba
 */

#include "NegfModel.h"
#include "Material.h"

NegfModel::NegfModel(const ModelOptions& options)
   : PhysicalModel(options)
{
}

NegfModel*
NegfModel::create(const Material* mat, const ModelOptions& options)
{
  return PhysicalModelInterface::create<NegfModel>(_create,_destroy, mat, options);
}

void
NegfModel::prepare_submodels(void)
{
    ModelOptions opts;
    opts.set_option("mass",1.0);
    create_submodels(_ham_models,"hamiltonian",opts);

}

