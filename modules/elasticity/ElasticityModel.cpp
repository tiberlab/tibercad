/*  
 * This file is part of the tiberCAD module elasticity.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file ElasticityModel.C
 * \brief tiberCAD elasticity module implementation.
 *
 * \note This file is part of module elasticity.
 */


#include "ElasticityModel.h"



using namespace std;


ElasticityModel*
ElasticityModel::create(const Material* mat, const ModelOptions& options)
{

  return PhysicalModel::create<ElasticityModel>(_create, _destroy, mat, options);

}


void
ElasticityModel::do_init(void)
{
  PhysicalModel::SubmodelIterator  it;
  it = submodels_begin("stiffness");
  StiffnessModel* _sm = dynamic_cast<StiffnessModel*> ((*it).second);
  _stiffness =  _sm->get_stiffness();


}

//! Print some useful information
void
ElasticityModel::do_print_info(void)
{
  PhysicalModel::SubmodelIterator it(submodels_begin("stiffness"));
  StiffnessModel* sm = dynamic_cast<StiffnessModel*> ((*it).second);
  sm->print_info();

}

void
ElasticityModel::calculate(const libMesh::Elem* elem, const libMesh::Point& point)
{
 
  //Update Body force model
  _force = 0;
  _strain = 0;
  _stress = 0;
  for (ID n = 0 ; n <_bfm.size() ; n++)
  {
    _bfm[n]->calculate(elem,point);
    _force  +=  _bfm[n]->get_force_source();
    _strain +=  _bfm[n]->get_strain_source();
    _stress +=  _bfm[n]->get_stress_source();
  }

}


void
ElasticityModel::prepare_submodels(void)
{
  ModelOptions opts;
  opts.set_option("type", "anisotropic");
  PhysicalModel* pm;
  create_submodel(pm, "stiffness", opts);

  create_submodels(_bfm, "body_force");
}
