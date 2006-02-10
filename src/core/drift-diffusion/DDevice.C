// $Id$

#include "DDevice.h"
#include "ElementData.h"
#include "SemiconductorModel.h"
#include "BoundaryData.h"
#include "BoundaryDescriptor.h"

// libmesh includes
#include "mesh.h"
#include "elem.h"



DD::Device::Device(Mesh* mesh, ElementData* element_data,
    BoundaryData* boundary_data)
{
  _mesh = mesh;
  _elem_data = element_data;
  _bound_data = boundary_data;

  // fill _materials for easier access to the different materials
  _materials = new MaterialList;
  ElementData::const_iterator it = _elem_data->begin();
  const ElementData::const_iterator end = _elem_data->end();
  for ( ; it != end; ++it)
    _materials->insert(it->second);

  // the same for boundary conditions
  _boundaries = new BoundaryList;
  if (_bound_data->get_size() > 0)
  {
    BoundaryData::const_iterator it = _bound_data->sides_begin();
    const BoundaryData::const_iterator end = _bound_data->sides_end();
    for ( ; it != end; ++it)
      _boundaries->insert(it->second);
  }
}

DD::Device::~Device(void)
{
  delete _materials;
  delete _boundaries;
}

bool
DD::Device::check_integrity(void) const
{
  bool is_ok = true;
  
  ElementData::const_iterator data_end = _elem_data->end();
  const_material_iterator material_end = _materials->end();
  
  MeshBase::const_element_iterator it =
    _mesh->level_elements_begin(0);
  const MeshBase::const_element_iterator end =
    _mesh->level_elements_end(0);

  if (it == end)
    is_ok = false;
  
  while (it != end)
  {
    ElementData::const_iterator el_it = _elem_data->find(*it);
    if (el_it == data_end)
    {
      is_ok = false;
      break;
    }
    if (_materials->find(el_it->second) == material_end)
    {
      is_ok = false;
      break;
    }
    ++it;
  }

  return is_ok;
}

BoundaryDescriptor*
DD::Device::get_boundary(const std::string& name) const
{
  BoundaryDescriptor* boundary = NULL;

  BoundaryData::const_iterator it = _bound_data->sides_begin();
  const BoundaryData::const_iterator end = _bound_data->sides_end();
  while ((it != end) && (it->second->get_id() != name))
    ++it;

  if (it != end)
    boundary = it->second;

  return boundary;
}


