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
 * \file MaterialBoundary.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/physics/MaterialBoundary.h"
#include "tibercad/physics/Material.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/io/Database.h"

#include <cassert>

using namespace std;

MaterialBoundary::MaterialBoundary(const ModelOptions& options)
 : PhysicalObject(BOUNDARY, options),
   _id_A(INVALID_ID),
   _id_B(INVALID_ID),
   _mat_A(NULL),
   _mat_B(NULL)
{
}


MaterialBoundary*
MaterialBoundary::create(ID id_A, Material* mat_A,
    ID id_B, Material* mat_B, const ModelOptions& options)
{
  // if it is an outer boundary, one of the two materials is NULL
  if (mat_A == NULL)
  {
    if (mat_B == NULL)
      throw InitFailedException("Trying to create boundary object touching "
          "no regions at all! Check your boundary regions.");

    std::swap(mat_A, mat_B);
    std::swap(id_A, id_B);
  }
  assert(mat_A != NULL);

  MaterialBoundary* mat = new MaterialBoundary(options);

  assert(mat != NULL);

  // If it is an internal boundary, then we create a new
  // database, otherwise we copy that of material A
  Database db;
  string name(mat_A->get_name());
  if (mat_B != NULL)
  {
    name += "%" +  mat_B->get_name();
    db.set_material(name, options.get_option("datafile", ""));
  }
  else
    db = mat_A->get_database();

  db.set_section("");
  mat->set_database(db);

  mat->_id_A = id_A;
  mat->_id_B = id_B;
  mat->_mat_A = mat_A;
  mat->_mat_B = mat_B;

  return mat;
}


void
MaterialBoundary::do_init(void)
{
  // we call the interface model initializer method
  ModelMap::iterator it(models_begin());
  ModelMap::const_iterator end(models_end());

  for ( ; it != end; ++it)
  {
    // the initialization might need a material
    (it->second)->set_material(_mat_A);
    (it->second)->init_interface(_mat_A, _mat_B);
  }
}
