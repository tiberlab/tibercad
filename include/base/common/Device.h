// $Id$

#ifndef __DEVICE_H__
#define __DEVICE_H__


#include "TypeDefs.h"

#include <vector>

class Material;
class Mesh;
class EquationSystems;

//! Higher-level definition of the  structure to  be  simulated.
/*!
 * This class contains all Material instances and the list of boundary
 * nodes as read from the mesh file
 */
class Device
{

  public:

    //! Constructor
    /*!
     * The mesh is assumed to be correctly prepared, i.e. all elements the
     * right subdomain id assigned (which is the physical region number
     * found in the mesh file).
     * The \c boundary_nodes map has to contain all nodes for each boundary
     * for which a boundary condition is implied.
     */
    Device(Mesh& mesh, BoundaryNodeMap& boundary_nodes);

    //! Destructor
    ~Device();

    //! Get a reference to the mesh
    Mesh& get_mesh(void) const;

    //! Get a reference to the equation systems object
    EquationSystems& get_equation_systems(void) const;

    //! Initialize this device
    /*!
     * This method will call \c init() of all Materials in this device
     */
    void init(void);

  
    //! Set a material for a geometrical region
    /*!
     * \c region_id is assumed to be a valid region number as given in the
     * mesh.
     *
     * \param material a pointer to the material
     * \param region_id the region number this material should belong to
     */
    void set_material(Material* material, ID region_id);
  
    //! Set a material for a number of geometrical regions
    /*!
     * \c region_id is assumed to be a vector of valid region numbers
     * as given in the mesh.
     *
     * \param material a pointer to the material
     * \param region_id the region numbers this material should belong to
     */
    void set_material(Material* material, const std::vector<ID>& region_ids);


    //! Get the material for a given region ID
    /*!
     * If \c region_id is not found in the material list, the NULL
     * pointer is returned.
     *
     * \param region_id the region number
     */
    const Material* get_material(ID region_id) const;
    
    //! Get the map that contains all boundary nodes for all boundaries
    BoundaryNodeMap& get_boundary_node_map(void) const;
    

  private:

    //! A typdef for convenience
    typedef std::map<ID, Material*> MaterialMap;

    //! The map that connects region number to material
    MaterialMap _material_map;

    //! The mesh for this device
    Mesh* _mesh;

    //! The equation systems used for this device
    /*!
     * This is stored here because it has to be consistent with the mesh
     * of the device
     */
    EquationSystems* _eq_system;

    //! A map that contains all nodes for boundary conditions
    BoundaryNodeMap* _boundary_nodes;
  
};


inline
const Material*
Device::get_material(ID region_id) const
{
  MaterialMap::const_iterator it(_material_map.find(region_id));

  if (it == _material_map.end())
    return NULL;

  return it->second;
}


inline
Mesh&
Device::get_mesh(void) const
{
  return *_mesh;
}


inline
EquationSystems&
Device::get_equation_systems(void) const
{
  return *_eq_system;
}


inline
BoundaryNodeMap&
Device::get_boundary_node_map(void) const
{
  return *_boundary_nodes;
}

#endif //  __DEVICE_H__
