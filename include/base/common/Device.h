// $Id$

#ifndef __DEVICE_H__
#define __DEVICE_H__


#include "TiberCad.h"
#include "TypeDefs.h"
#include "ElementSide.h"
#include "HashMap.h"
#include "IDSet.h"
#include "ModelOptions.h"
#include "DeviceException.h"

#include "elem.h"

#include <vector>
#include <set>

class Material;
class MaterialBoundary;
class EdgeObject;
class NodeObject;
class MeshBase;
class EquationSystems;
class AtomisticStructure;
class MeshRegionInfo;
class BoundaryRegions;
class Point;
class QuantumContact;

//! Higher-level definition of the  structure to  be  simulated.
/*!
 * This class contains all Material instances and the list of boundary
 * nodes as read from the mesh file
 */
class Device
{

  private:

    //! A typedef for the map containing atomistic structures
    typedef std::map<std::string,  AtomisticStructure*> AtomStructMap;

    typedef std::map<std::string, QuantumContact*> QuantumContactMap;


  public:

    typedef std::map<ID, std::vector<ID> > BCNodeMap;

    //! The atomistic structures iterator type
    typedef AtomStructMap::iterator atomistic_structure_iterator;

    //! The quantum contact iterator type
    typedef QuantumContactMap::iterator quantum_contact_iterator;


    //! Destructor
    ~Device();

    //! The method for creation of a device
    /*!
     * \param options the options needed for device creation
     * \return a pointer to the newly created device
     *
     * \c options has to contain at least the following options:
     * \li "meshfile" -> filename
     * \li "mesh_units" -> the units of the mesh (cf. get_mesh_units())
     */
    static Device* create(const ModelOptions& options) TBDLLOCAL;


    //! Destroy a Device object
    static void destroy(Device* device) TBDLLOCAL;


    //! Get a reference to the mesh
    MeshBase& get_mesh(void) const;


    //! Get a reference to the equation systems object
    EquationSystems& get_equation_systems(void) const;


    //! Get a reference to the equation systems object for a given mesh
    /*!
     * If it does not exist yet, it will be created.
     */
    EquationSystems& get_equation_systems(MeshBase* mesh);


    //! Prepare the device
    /*!
     * This creates all device regions, materials and
     * atomistic structures.
     */
    void prepare(void);


    //! Initialize this device
    /*!
     * This method will call \c init() of all Materials in this device
     */
    void init(void) TBDLLOCAL;


    //! Set a material for a number of geometrical regions
    /*!
     * \c region_id is assumed to be a vector of valid region numbers
     * as given in the mesh.
     *
     * \param material a pointer to the material
     * \param region_id the region numbers this material should belong to
     * \param region_name the name to be assigned
     */
    void set_material(Material* material, const std::vector<ID>& region_ids,
        const std::string& region_name) TBDLLOCAL;


    //! Get the material for a given region ID
    /*!
     * If \c region_id is not found in the material list, the NULL
     * pointer is returned.
     *
     * \param region_id the region number
     */
    const Material* get_material(ID region_id) const;


    /*! \copydoc get_material(ID) const */
    Material* get_material(ID region_id);


    //! Get the material for a region name
    /*!
     * \param name the user defined name of a region
     * \return the material pointer or \c NULL if \c name does not name
     * a valid region
     */
    const Material* get_material(const std::string& name) const;


    /*! \copydoc get_material(const std::string&) const */
    Material* get_material(const std::string& name);


    //! Get the material for a given element
    const Material* get_material(const Elem* elem) const;


    /*! \copydoc get_material(const Elem*) const */
    Material* get_material(const Elem* elem);


    //! Get the MaterialBoundary object for a given ID
    /*!
     * if there is no object associated to \c id
     * it will be created.
     */
    MaterialBoundary* get_boundary_object(ID id);


    //! Get the MaterialBoundary object for a given element side
    MaterialBoundary* get_boundary_object(const Elem*, int side);


    //! Get the EdgeObject for a given ID
    /*!
     * if there is no object associated to \c id
     * it will be created.
     */
    EdgeObject* get_edge_object(ID id);


    //! Get the EdgeObject for a given element edge
    EdgeObject* get_edge_object(const Elem*, int edge);


    //! Get the NodeObject for a given ID
    /*!
     * if there is no object associated to\c id
     * it will be created.
     */
    NodeObject* get_node_object(ID id);


    //! Get the NodeObject for a given element node
    NodeObject* get_node_object(const Elem*, int node);


    //! Get the NodeObject for a given node
    NodeObject* get_node_object(const Node* node);


    //! Get the map that contains all boundary nodes for all boundaries
    BCNodeMap& get_boundary_node_map(void) const;


    //! Get the mesh units
    /*!
     * Mesh units are in SI units, i.e. meters. The return value is
     * the distance in real space corresponding to a distance of 1 in
     * the mesh object, or:
     * \f[\mathrm{d}_{real}(x, y) = \gamma\mathrm{d}_{mesh}(x,y)\f]
     * where \f$\gamma\f$ is the result of get_mesh_units()
     */
    double get_mesh_units(void) const;


    //! Get the set with all region IDs
    /*!
     * Only active regions are returned (that is regions that have a material
     * associated).
     */
    const IDSet& get_active_region_ids(void) const TBDLLOCAL;


    //! Get the name of a region
    /*!
     * \param id the ID of the physical region
     *
     * If a region has no name associated, it will be assigned the
     * empty string.
     */
    const std::string& get_region_name(ID id) const;


    // Get the name of a boundary region
    /*
     * \param id the ID of the boundary region
     *
     * If a region has no name associated, it will be assigned the
     * empty string.
     */
    //const std::string& get_boundary_region_name(ID id) const;


    //! Get the region IDs of the region with name \c name
    /*!
     * \c name can be the name of a physical region, of a cluster
     * or of a mesh region
     * Only active regions are returned (that is regions that have a material
     * associated).
     */
    void get_active_region_ids(const std::string& name, std::set<ID>& ids) const TBDLLOCAL;


    //! Get the region IDs of the mesh region with name \c name
    /*!
     * This looks only in the list of original mesh region names.
     */
    void get_mesh_region_ids(const std::string& name, std::vector<ID>& ids) const TBDLLOCAL;


    //! Extract physical regions from a string
    void extract_physical_regions(const std::string& str, IDSet& ids) const;


    /* // ! Get a const reference to the boundary region descriptor */
    //const BoundaryRegions& get_boundary_regions(void) const;


    //! Get the region IDs of the boundary region with name \c name
    /*!
     * \c name can be one of
     * \li boundary region name
     * \li material interface specification of type \c matA/matB, e.g. \c Si/Ge
     * \li region interface specification of type \c regA/regB
     * \li a mix of the latter two
     *
     * \note An interface specification with the two components being the same
     * is ignored.
     */
    void get_boundary_region_ids(const std::string& name,
        std::vector<ID>& ids) const;


    //! Set the name for a boundary region
    //void set_boundary_region_name(const std::string& name,
    //    const std::vector<ID>& ids);


    //! Define a cluster
    void set_cluster(const std::string& name, const std::vector<ID>& ids) TBDLLOCAL;


    //! Get the type of symmetry
    TiberCad::Symmetry get_symmetry(void) const;


    //! Get the pointer for the atomistic structure defined with \c name
    AtomisticStructure* get_atomistic_structure(const std::string& name);

    //! Get the requested atomistic structures in a vector
    void get_atomistic_structures(const std::string& names,
        std::vector<AtomisticStructure*>& structures);


    //! Get the iterator to the first atomistic structure
    atomistic_structure_iterator atomistic_structures_begin(void);

    //! Get the past-the-end iterator for the atomistic structures
    atomistic_structure_iterator atomistic_structures_end(void);


    QuantumContact* get_quantum_contact(const std::string& name);

    //! Get the iterator to the first quantum contact
    quantum_contact_iterator quantum_contacts_begin(void);

    //! Get the past-the-end iterator for the quantum contacts
    quantum_contact_iterator quantum_contacts_end(void);



  private:


    //! A typdef for the bulk materials
    typedef HashMap<ID, Material*>::Type MaterialMap;

    //! A typdef for the material boundaries
    typedef HashMap<ID, MaterialBoundary*>::Type BoundaryMap;

    //! A typdef for the edge objects
    typedef HashMap<ID, EdgeObject*>::Type EdgeObjMap;

    //! A typdef for the node objects
    typedef HashMap<ID, NodeObject*>::Type NodeObjMap;


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
    Device(void) TBDLLOCAL;


    //! Set a material for a physical region
    /*!
     * \param material a pointer to the material
     * \param region_id the region number this material should belong to
     * \throws {InitFailedException if \c region_id is invalid or already
     * used.}
     */
    void set_material(Material* material, ID region_id) TBDLLOCAL;


    //! Set the name for a physical region
    void set_region_name(const std::string& name, const std::vector<ID>& ids) TBDLLOCAL;


    //! Set options for this device
    /*!
     * The options are stored internally and are accessible through
     * special methods.
     * Options have to be specified at creation time.
     */
    void set_options(const ModelOptions& options) TBDLLOCAL;


    //! Get a reference to the model options
    ModelOptions& get_options(void) TBDLLOCAL;


    //! Creates the mesh and the equation system
    /*!
     * This method assumes, that \c _options contain the name of the
     * meshfile and the dimension
     */
    void setup_mesh(void) TBDLLOCAL;


    //! Setup all regions and materials
    void setup_regions(void);


    //! Setup clusters
    /*!
     * A region cluster contains different physical regions of the mesh
     * which possibly overlap with the material regions described in the
     * \c Region sections
     */
    void setup_clusters(void);


    //! Setup atomistic structures
    void setup_atomistic_structures(void);

    //! Setup quantum contacts
    void setup_quantum_contacts(void);

    //! Prepares the boundaries
    /*!
     * Assign to each explicitly or implicitly (region interfaces) defined
     * boundary a unique index.
     */
    //void prepare_boundaries(void) TBDLLOCAL;


    //! The map that connects region number to material
    MaterialMap _material_map;


    //! The map connecting boundary regions to model containers
    BoundaryMap _boundary_map;


    //! The map connecting edge regions to model containers
    EdgeObjMap _edge_map;


    //! The map connecting node regions to model containers
    NodeObjMap _node_map;



    //! The map that connects atomistic structure names to pointers
    /*! (keep track of existing atomistic structures) */
    AtomStructMap _atomistic_structure_map;


    //! The mesh for this device
    MeshBase* _mesh;


    //! A typdef for the bulk materials
    //typedef HashMap<Elem*, std::vector<unsigned int>* >::Type ElemAtomsMap;


    //! A map between elements of _mesh and atoms
    //! It is defined as vector for fast indexing on elem ID.
    // std::map<AtomisticStructure*, ElemAtomsMap> _as_elematom_map;


    //! The mesh unit in m
    /*!
     * A distance of 1 in the mesh corresponds to \c _mesh_units m
     */
    double _mesh_units;


    //! The equation systems used for this device
    /*!
     * This is stored here because it has to be consistent with the mesh
     * of the device
     */
    EquationSystems* _eq_system;


    //! A map with equations systems for all used meshes
    std::map<const MeshBase*, EquationSystems*> _eq_sys_map;


    //! A map that contains all nodes for boundary conditions
    BCNodeMap* _boundary_nodes;


    //! User defined options for this device
    ModelOptions _options;


    //! A set with all region IDs
    //std::set<ID> _region_ids;


    //! A set with all active region IDs
    /*!
     * An active region is a region with an associated material.
     */
    IDSet _active_region_ids;


    //! A map that assigns physical region or cluster IDs to names
    std::map<ID, std::string> _region_names;


    //! A structure containing the original mesh region info
    MeshRegionInfo* _mesh_region_info;

    //! A structure containing the original boundary region info
    BoundaryRegions* _bd_regions;

    //! The map that connects quantum contacts to pointers
     /*! (keep track of existing quantum contacts) */
    QuantumContactMap _quantum_contact_map;

    //! A map that assigns boundary region IDs to boundary region names
    //std::map<ID, std::string> _boundary_region_names;


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
const Material*
Device::get_material(const Elem* elem) const
{
  return get_material(elem->subdomain_id());
}


inline
Material*
Device::get_material(const Elem* elem)
{
  return get_material(elem->subdomain_id());
}


inline
MeshBase&
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
Device::BCNodeMap&
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
const IDSet&
Device::get_active_region_ids(void) const
{
  return _active_region_ids;
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


/*
inline
const BoundaryRegions&
Device::get_boundary_regions(void) const
{
  return *_bd_regions;
}
*/

/*
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
*/


inline
TiberCad::Symmetry
Device::get_symmetry(void) const
{
  return _symmetry;
}


inline
Device::atomistic_structure_iterator
Device::atomistic_structures_begin(void)
{
  return _atomistic_structure_map.begin();
}

inline
Device::atomistic_structure_iterator
Device::atomistic_structures_end(void)
{
  return _atomistic_structure_map.end();
}


inline
Device::quantum_contact_iterator
Device::quantum_contacts_begin(void)
{
  return _quantum_contact_map.begin();
}

inline
Device::quantum_contact_iterator
Device::quantum_contacts_end(void)
{
  return _quantum_contact_map.end();
}



inline
AtomisticStructure*
Device::get_atomistic_structure(const std::string& name)
{
  if (_atomistic_structure_map.find(name) != _atomistic_structure_map.end())
    return _atomistic_structure_map[name];
  else
    return NULL;
}

inline
QuantumContact*
Device::get_quantum_contact(const std::string& name)
{
  if (_quantum_contact_map.find(name) != _quantum_contact_map.end())
    return _quantum_contact_map[name];
  else
    return NULL;
}

#endif //  __DEVICE_H__
