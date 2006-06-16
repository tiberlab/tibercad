// $Id$

#ifndef _DDEVICE_H_
#define _DDEVICE_H_

#include <vector>
#include <set>
#include <map>


// forward declarations
class Mesh;
class ElementData;
class BoundaryData;
class DriftDiffusionProperties;
class ElectricalContact;

namespace DD {
class Device
{
  public:

    typedef std::set<DriftDiffusionProperties*> MaterialList;
    typedef MaterialList::iterator material_iterator;
    typedef MaterialList::const_iterator const_material_iterator;

    typedef std::set<ElectricalContact*> BoundaryList;
    typedef BoundaryList::iterator boundary_iterator;
    typedef BoundaryList::const_iterator const_boundary_iterator;

    /**
     * Constructor
     * 
     * The constructor needs a valid mesh and consistent \p ElementData
     * and \p BoundaryData structures.
     * Boundary data and Element data have to be given for level 0
     * elements. Any operation on the mesh therefore has to maintain
     * level 0 structure.
     */
    Device(Mesh* mesh, ElementData* element_data,
        BoundaryData* boundary_data);

    ~Device(void);

    Mesh& get_mesh(void) const;

    ElementData& get_element_data(void) const;

    BoundaryData& get_boundary_data(void) const;

    MaterialList& get_materials(void) const;

    BoundaryList& get_boundaries(void) const;

    /**
     * @returns a pointer to the boundary descriptor of the given
     * boundary with name @c name.
     *
     * If @c name does not exist, @c NULL is returned
     */
    ElectricalContact* get_boundary(const std::string& name) const;

    /**
     * Checks integrity of \p this object
     *
     * Integrity means: every 0-level element of the mesh is
     * connected with a material.
     */
    bool check_integrity(void) const;


  private:

    Mesh* _mesh;
    ElementData*  _elem_data;
    BoundaryData* _bound_data;
    
    MaterialList*      _materials;
    BoundaryList*      _boundaries;


};


//
// inline member functions
// 

inline
Mesh&
Device::get_mesh(void) const
{
  return *_mesh;
}

inline
ElementData&
Device::get_element_data(void) const
{
  return *_elem_data;
}

inline
BoundaryData&
Device::get_boundary_data(void) const
{
  return *_bound_data;
}

inline
Device::MaterialList&
Device::get_materials(void) const
{
  return *_materials;
}

inline
Device::BoundaryList&
Device::get_boundaries(void) const
{
  return *_boundaries;
}

};

#endif // _DDEVICE_H_
