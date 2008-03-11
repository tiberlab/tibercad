// $Id$


#include "MeshUtils.h"

#include "mesh_data_elements.h"
#include "mesh.h"
#include "elem.h"


void
MeshUtils::assign_subdomain_ids(Mesh& mesh, MeshData_elements& meshdata)
{
  Mesh::element_iterator it = mesh.local_elements_begin();
  const Mesh::element_iterator end = mesh.local_elements_end();

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;

    if (meshdata.has_data(elem))
    {
      int id = static_cast<int>(meshdata(elem));
      elem->subdomain_id() = id;
    }
  }
}



void
MeshUtils::get_subdomain_ids(Mesh& mesh, std::set<ID>& subdomain_ids)
{
  subdomain_ids.clear();

  Mesh::element_iterator it = mesh.local_elements_begin();
  const Mesh::element_iterator end = mesh.local_elements_end();

  for ( ; it != end; ++it)
  {
    Elem* elem = *it;

    ID id = static_cast<ID>(elem->subdomain_id());
    subdomain_ids.insert(id);
  }
}
