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
  assert(mat_A != NULL);

  // TODO how to use database correctly in this case?

  string name(mat_A->get_name());
  if (mat_B != NULL)
    name += "%" +  mat_B->get_name();

  Database db(name, options.get_option("datafile", ""));
  db.set_section("");

  MaterialBoundary* mat = new MaterialBoundary(options);

  assert(mat != NULL);

  mat->_id_A = id_A;
  mat->_id_B = id_B;
  mat->_mat_A = mat_A;
  mat->_mat_B = mat_B;

  mat->set_database(db);

  return mat;
}


void
MaterialBoundary::do_init(void)
{
  // we call the interface model initializer method
  ModelMap::iterator it(models_begin());
  ModelMap::const_iterator end(models_end());

  for ( ; it != end; ++it)
    (it->second)->init_interface(_mat_A, _mat_B);
}
