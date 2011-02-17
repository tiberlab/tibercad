// $Id$

#include "MaterialBoundary.h"
#include "Material.h"
#include "PhysicalModel.h"
#include "Database.h"

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
