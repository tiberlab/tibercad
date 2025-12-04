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
 * \file AugmentSparsityForMassConservation.C
 * \brief tiberCAD dssc module implementation.
 *
 * \note This file is part of module dssc.
 */


#include "AugmentSparsityForMassConservation.h"

#include "libmesh/implicit_system.h"
#include "libmesh/node.h"

using namespace std;
using namespace libMesh;

AugmentSparsityForMassConservation::AugmentSparsityForMassConservation(
    ImplicitSystem& sys,
    const std::map<dof_id_type, std::set<unsigned int>>& coupled_vars)
:
  _sys(sys),
  _coupled_vars(coupled_vars)
{
}


void
AugmentSparsityForMassConservation::augment_sparsity_pattern(
    libMesh::SparsityPattern::Graph & ,
    std::vector<dof_id_type> & n_nz,
    std::vector<dof_id_type> & n_oz)
{
  //cerr << "=========> Augment sparsity pattern\n";

  const MeshBase& mesh = _sys.get_mesh();
  const DofMap& dof_map = _sys.get_dof_map();

  dof_id_type n_local_dofs = dof_map.n_local_dofs();
  dof_id_type n_remote_dofs = dof_map.n_dofs() - n_local_dofs;

  map<unsigned int, set<unsigned int>>::iterator it(_coupled_vars.begin());
  const map<unsigned int, set<unsigned int>>::iterator end(_coupled_vars.end());

  for ( ; it != end; ++it)
  {
    dof_id_type dof = it->first;
    //cerr << "    DOF id : " << dof << endl;

    if ((dof >= dof_map.first_dof()) &&
        (dof < dof_map.end_dof()))
    {
      n_nz[dof - dof_map.first_dof()] = n_local_dofs;
      n_oz[dof - dof_map.first_dof()] = n_remote_dofs;
    }
  }



//  MeshBase::const_node_iterator nd = mesh.local_nodes_begin();
//  MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();

//  for ( ; nd != nd_end; ++nd)
//  {

//  }
}
