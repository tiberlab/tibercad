/*  
 * This file is part of the tiberCAD module dssc.
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
 * \file AugmentSparsityForMassConservation.h
 * \brief tiberCAD dssc module header.
 *
 * \note This file is part of module dssc.
 */


#ifndef AUGMENTSPARSITYFORMASSCONSERVATION_H_
#define AUGMENTSPARSITYFORMASSCONSERVATION_H_


#include "libmesh/dof_map.h"

namespace libMesh {
  class ImplicitSystem;
}

using libMesh::DofMap;
using libMesh::ImplicitSystem;
using libMesh::dof_id_type;

/*!
 * \brief Augment the sparsity pattern for mass conservation
 */
class AugmentSparsityForMassConservation : public DofMap::AugmentSparsityPattern
{

  public:

    /*!
     * constructor
     */
    AugmentSparsityForMassConservation(ImplicitSystem& sys,
        const std::map<dof_id_type, std::set<unsigned int>>& coupled_vars);


    /*!
     * Our function to do the actual augment
     */
    virtual void augment_sparsity_pattern(libMesh::SparsityPattern::Graph & ,
                                          std::vector<dof_id_type> & n_nz,
                                          std::vector<dof_id_type> & n_oz);


  private:

    /*!
     * The equation systems object
     */
    ImplicitSystem& _sys;


    /*!
     * The other variables the DOF should couple to
     */
    std::map<unsigned int, std::set<unsigned int>> _coupled_vars;
};

#endif /* AUGMENTSPARSITYFORMASSCONSERVATION_H_ */
