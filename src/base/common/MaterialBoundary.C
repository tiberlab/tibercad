// $Id$

#include "MaterialBoundary.h"
#include "Material.h"

#include <cassert>

using namespace std;

MaterialBoundary::MaterialBoundary(Material* mat_A, Material* mat_B
    , const ModelOptions& options)
 : PhysicalObject(BOUNDARY, options),
   _mat_A(mat_A),
   _mat_B(mat_B)
{
}


MaterialBoundary*
MaterialBoundary::create(Material* mat_A, Material* mat_B,
    const ModelOptions& options)
{
  assert(mat_A != NULL);
  //string name(mat_A->get_name());
  //if (mat_B != NULL)
  //  name += "-" +  mat_B->get_name();

  MaterialBoundary* mat = new MaterialBoundary(mat_A, mat_B, options);

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
    // mat_A is never NULL
    PhysicalModel* pma = _mat_A->get_model(it->first);
    PhysicalModel* pmb = NULL;
    if (_mat_B != NULL)
      pmb = _mat_B->get_model(it->first);

    //(it->second)->init_boundary(pma, pmb);
  }
}
