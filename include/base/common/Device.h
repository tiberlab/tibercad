// $Id$

#ifndef __DEVICE_H__
#define __DEVICE_H__


#include "TiberCad.h"
#include "TypeDefs.h"
#include "ModelOptions.h"
#include "DeviceException.h"

#include "elem.h"
#include "mesh_data_elements.h"

#include <vector>
#include <set>

class Material;
class Control;
class Mesh;
class EquationSystems;
class AtomisticStructure;

//! Higher-level definition of the  structure to  be  simulated.
/*!
 * This class contains all Material instances and the list of boundary
 * nodes as read from the mesh file
 */
class Device
{

  public:

    //! Destructor
    ~Device();

    //! The method for creation of a device
    /*!
     * \param options the options needed for device creation
     * \return a pointer to the newly created device
     *
     * \c options has to contain at least the following options:
     * \li "meshfile" -> filename
     * \li "dimension" -> the real space dimension (1, 2 or 3)
     * \li "mesh_units" -> the units of the mesh (cf. get_mesh_units())
     */
    static Device* create(const ModelOptions& options);


    //! Destroy a Device object
    static void destroy(Device* device);


    //! Get a reference to the mesh
    Mesh& get_mesh(void) const;

    //! Get a pointer to the meshdata
    MeshData_elements* get_meshdata(void) const;
    
    //! Get a reference to the equation systems object
    EquationSystems& get_equation_systems(void) const;


    //! Initialize this device
    /*!
     * This method will call \c init() of all Materials in this device
     */
    void init(void);

  
    //! Set a material for a physical region
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


    //! Set an atomistic structure to be kept in structures map
    void set_atomistic_structure(const std::string& name,
        AtomisticStructure* atomistic_structure);


    //! Set the Control module which will control this device
    void set_control(Control* control);


    /*! \copydoc get_material(ID) const */
    Material* get_material(ID region_id);

    
    //! Get the material for a given region ID
    /*!
     * If \c region_id is not found in the material list, the NULL
     * pointer is returned.
     *
     * \param region_id the region number
     */
    const Material* get_material(ID region_id) const;


    /*! \copydoc get_material(const std::string&) const */
    Material* get_material(const std::string& name);


    //! Get the material for a region name
    /*!
     * \param name the user defined name of a region
     * \return the material pointer or \c NULL if \c name does not name
     * a valid region
     */
    const Material* get_material(const std::string& name) const;
    
    
    //! Get the map that contains all boundary nodes for all boundaries
    BoundaryNodeMap& get_boundary_node_map(void) const;

    
    //! Get the mesh units
    /*!
     * Mesh units are in SI units, i.e. meters. The return value is
     * the distance in real space corresponding to a distance of 1 in
     * the mesh object, or:
     * \f[\mathrm{d}_{real}(x, y) = \gamma\mathrm{d}_{mesh}(x,y)\f]
     * where \f$\gamma\f$ is the result of get_mesh_units()
     */
    double get_mesh_units(void) const;


    //! Get a reference to the Control module that controls this device
    Control& get_control(void) const;


    //! Get the set with all region IDs
    const std::set<ID>& get_region_ids(void) const;


    //! Get the name of a region
    /*!
     * \param id the ID of the physical region
     *
     * If a region has no name associated, it will be assigned the
     * empty string.
     */
    const std::string& get_region_name(ID id) const;


    //! Get the name of a boundary region
    /*!
     * \param id the ID of the boundary region
     *
     * If a region has no name associated, it will be assigned the
     * empty string.
     */
    const std::string& get_boundary_region_name(ID id) const;


    //! Get the region IDs of the region with name \c name
    /*!
     * \c name can be the name of a region or of a cluster
     */
    void get_region_ids(const std::string& name, std::vector<ID>& ids) const;


    //! Get the region IDs of the boundary region with name \c name
    void get_boundary_region_ids(const std::string& name,
        std::vector<ID>& ids) const;


    //! Set the name for a region
    void set_region_name(const std::string& name, const std::vector<ID>& ids);


    //! Set the name for a boundary region
    void set_boundary_region_name(const std::string& name,
        const std::vector<ID>& ids);


    //! Define a cluster
    void set_cluster(const std::string& name, const std::vector<ID>& ids);


    //! Get the type of symmetry
    TiberCad::Symmetry get_symmetry(void) const;
    

    //! Get the pointer for the atomistic structure defined with \c name
    AtomisticStructure* get_atomistic_structure(const std::string&);


  private:


    //! A typdef for convenience
    typedef std::map<ID, Material*> MaterialMap;

    //! A typdef for convenience
    typedef std::map<std::string, std::vector<ID> > ClusterMap;


    //! Empty Constructor
    /*!
     * The mesh is assumed to be correctly prepared, i.e. all elements the
     * right subdomain id assigned (which is the physical region number
     * found in the mesh file).
     * The \c boundary_nodes map has to contain all nodes for each boundary
     * for which a boundary condition is implied.
     */
    Device(void);

    
    //! Set options for this device
    /*!
     * The options are stored internally and are accessible through
     * special methods.
     * Options have to be specified at creation time.
     */
    void set_options(const ModelOptions& options);

    
    //! Get a reference to the model options
    ModelOptions& get_options(void);

    
    //! creates the mesh and the equation system
    /*!
     * This method assumes, that \c _options contain the name of the
     * meshfile and the dimension
     */
    void setup_mesh(void);

    
    //! The map that connects region number to material
    MaterialMap _material_map;


    //! The map that connects atomistic structure names to pointers
    /*! (keep track of existing atomistic struxctures) */
    std::map<std::string,  AtomisticStructure* > _atomistic_structure_map;
  

    //! The mesh for this device
    Mesh* _mesh;

    //! The meshdata for this device
    MeshData_elements* _meshdata;

   
    //! The mesh unit in m
    /*!
     * A distance of 1 in the mesh corresponds to \c _mesh_units m
     */
    double _mesh_units;


    //! The control module which controls this device
    Control* _control;

    
    //! The equation systems used for this device
    /*!
     * This is stored here because it has to be consistent with the mesh
     * of the device
     */
    EquationSystems* _eq_system;

    
    //! A map that contains all nodes for boundary conditions
    BoundaryNodeMap* _boundary_nodes;


    //! User defined options for this device
    ModelOptions _options;


    //! A set with all region IDs
    std::set<ID> _region_ids;
  

    //! A map that assigns region IDs to region names
    std::map<ID, std::string> _region_names;
  

    //! A map that assigns boundary region IDs to boundary region names
    std::map<ID, std::string> _boundary_region_names;


    //! A map containing all clusters
    ClusterMap _cluster_map;


    //! The symmetry of the device
    /*!
     * The default assumes no special symmetry.
     */
    TiberCad::Symmetry _symmetry;

};





//
// inline methods
// 

inline
void
Device::destroy(Device* device)
{
  delete device;
}


inline
void
Device::set_options(const ModelOptions& options)
{
  _options = options;
}


inline
ModelOptions&
Device::get_options(void)
{
  return _options;
}


inline
Material*
Device::get_material(ID region_id)
{
  MaterialMap::const_iterator it(_material_map.find(region_id));

  if (it == _material_map.end())
    return NULL;

  return it->second;
}


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
MeshData_elements*
Device::get_meshdata(void) const
{
  return  _meshdata;
}



inline
void
Device::set_control(Control* control)
{
  _control = control;
}


inline
Control&
Device::get_control(void) const
{
  return *_control;
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

inline
double
Device::get_mesh_units(void) const
{
  return _mesh_units;
}


inline
const std::set<ID>&
Device::get_region_ids(void) const
{
  return _region_ids;
}



inline
const std::string&
Device::get_region_name(ID id) const
{
  std::map<ID, std::string>::const_iterator it(_region_names.find(id));

  if (it == _region_names.end())
  {
    std::ostringstream s;
    s << "Tried to access unknown region with id " << id;
    throw (DeviceException(s.str()));
  }
  return it->second;
}




inline
const std::string&
Device::get_boundary_region_name(ID id) const
{
  std::map<ID, std::string>::const_iterator it(_boundary_region_names.find(id));

  if (it == _boundary_region_names.end())
  {
    std::ostringstream s;
    s << "Tried to access unknown boundary region with id " << id;
    throw (DeviceException(s.str()));
  }
  return it->second;
}



inline
TiberCad::Symmetry
Device::get_symmetry(void) const
{
  return _symmetry;
}



inline
void 
Device::set_atomistic_structure(const std::string& name, AtomisticStructure* atomistic_structure)
{
  _atomistic_structure_map[name] = atomistic_structure;
}


inline
AtomisticStructure*
Device::get_atomistic_structure(const std::string& name)
{
  return _atomistic_structure_map[name];
} 
    
#endif //  __DEVICE_H__
