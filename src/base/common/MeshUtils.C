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


bool MeshUtils::may_belong_to_element(const Elem* element, Point& point)
{
  const unsigned int n = element->n_nodes();
  double  min_x;double  min_y; double  min_z;
  double  max_x;double  max_y; double  max_z;

  Point vertex = element->point(0);
  min_x = vertex(0); min_y = vertex(1); min_z = vertex(2);
  max_x = min_x;  max_y = min_y; max_z = min_z;

  for (unsigned int i = 1 ; i < n ; i++)
    {
      vertex = element->point(i);
      double x = vertex(0);
      double y = vertex(1);
      double z = vertex(2);

      if (min_x > x) min_x = x; if (min_y > y) min_y = y; if (min_z > z) min_z = z;

      if (max_x < x) max_x = x; if (max_y < y) max_y = y; if (max_z < z) max_z = z;



    }

  if ( (point(0) > max_x) ||  (point(0) < min_x) ||
       (point(1) > max_y) ||  (point(1) < min_y) ||
       (point(2) > max_z) ||  (point(2) < min_z) )

    {    return(false);}

  else
    {    return(true) ; }



}

