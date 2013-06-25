// $Id$

#ifndef _SIMULATIONENVIRONMENT_H_
#define _SIMULATIONENVIRONMENT_H_

#include "TypeDefs.h"
#include "HashSet.h"
#include "Device.h"
#include "BoundaryNodeMap.h"
#include "BoundaryElementMap.h"
#include "tiber_dll.h"

#include "elem.h"

#include <set>
#include <cassert>

class Boundary;
class Node;

//! Contains special data needed by a solver
/*!
 * It contains a reference to the device, a list of elements which define
 * the region for this simulation, a list of boundary nodes, the boundary
 * conditions.
 */
class SimulationEnvironment
{

  private:

    //! A typedef for convenience
    typedef std::set<Boundary*> BCMap;

    //! A typedef for convenience
    typedef HashMap<ElementSide, ID, ElementSide::hash>::Type ElemSideMap;

    //! The type of the element list
    typedef HashSet<const Elem*>::Type ElementList;


  public:

    //! An iterator for the element sides lying on a boundary
    class BoundarySideIterator;

    //! An iterator for the nodes lying on a boundary
    typedef std::map<const Node*, ID>::const_iterator BoundaryNodeIterator;

    //! An iterator for the elements
    typedef ElementList::const_iterator ConstElemIterator;

    //! An iterator for the region IDs
    typedef std::set<ID>::const_iterator RegionIDIterator;

    //! An iterator for the boundaries
    typedef std::set<Boundary*>::iterator BoundaryIterator;
    /*
    class BoundaryIterator
    {
      public:

        BoundaryIterator(void) {};

        BoundaryIterator(const BCMap& bd, const BCMap::const_iterator it) :
          _bd_it(it), _bd_end(bd.end()) {}

        BoundaryIterator(const BoundaryIterator& it) :
          _bd_it(it._bd_it), _bd_end(it._bd_end) {}

        BoundaryIterator& operator++(void)
        {
          Boundary* bd = _bd_it->second;
          while ((_bd_it != _bd_end) && (bd == _bd_it->second))
            ++_bd_it;

          return *this;
        }

        BoundaryIterator& operator=(const BoundaryIterator& it)
        {
          _bd_it = it._bd_it;
          _bd_end = it._bd_end;
        }

        bool operator==(const BoundaryIterator& it)
        {
          return ((_bd_it == it._bd_it) || (_bd_it->second == it._bd_it->second));
        }

        bool operator!=(const BoundaryIterator& it)
        {
          return !(*this == it);
        }

        Boundary* operator*(void)
        {
          return ((_bd_it == _bd_end) ? NULL : _bd_it->second);
        }


      private:

        BCMap::const_iterator _bd_it;
        BCMap::const_iterator _bd_end;

    };
    */

    //! The constructor
    /*!
     * \param device the device
     * \param region_numbers the regions to be used for this simulation
     */
    SimulationEnvironment(Device& device,
        const std::set<ID>& region_numbers);


    //! Destructor
    ~SimulationEnvironment(void);


    //! Get a reference to the device
    Device& get_device(void);


    //! Get a const reference to the device
    const Device& get_device(void) const;


    //! Get a reference to the mesh
    MeshBase& get_mesh(void);


    //! Set the mesh
    void set_mesh(MeshBase* mesh);



    //! Add a boundary
    /*!
     * \param boundary the boundary to assign to this
     * simulation
     */
    void add_boundary(Boundary* boundary) TBDLLOCAL;


    //! Prepare structures that are needed for other setup
    /*!
     * Creates a list with all elements belonging to this simulation
     */
    void prepare(void) TBDLLOCAL;


    //! Initialize the environment
    /*!
     * Initialization does the following things:
     * \li create data structures to get access to the boundary properties
     */
    void init(void) TBDLLOCAL;


    //! Prepares the environment for a solve
    /*!
     * This sets all elements for this simulation to active and the others
     * to inactive.
     */
    void prepare_for_solve(void);


    //! Check if a node is part of a given boundary
    bool is_node_on_boundary(const Node* node, const Boundary* boundary) const;


    //! Get the boundary for a given boundary number
    /*!
     * \param boundary_number the number of the boundary
     * \return a pointer to the boundary if found, \c NULL otherwise
     */
    Boundary* get_boundary(ID boundary_number) const;


    //! Get the boundary for a given ElementSide
    /*!
     * \param side an element side
     * \return a pointer to the boundary if found, \c NULL otherwise
     */
    Boundary* get_boundary(const ElementSide& side) const;


    //! Get the boundary with user defined name \c name
    /*!
     * \param name the user defined name to look for
     * \return a pointer to the boundary if found, \c NULL otherwise
     */
    Boundary* get_boundary(const std::string& name) const;


    //! Get the boundary for a given Node
    /*!
     * \param node a node
     * \return a pointer to the boundary if found, \c NULL otherwise
     *
     * \note If two boundaries touch in node \c node, the return value
     * is unpredictable.
     */
    Boundary* get_boundary(const Node* node) const;


    //! Get the list of all Nodes that belong to a certain boundary
    /*!
     * \param name the name of the boundary to look for
     * \param nodelist the list to create
     */
    void get_boundary_nodes(const std::string& name,
        std::set<const Node*>& nodelist);


    //! Get the list of all Nodes that belong to a certain boundary
    /*!
     * \param boundary the pointer to the boundary to look for
     * \param nodelist the list to create
     */
    void get_boundary_nodes(const Boundary* boundary,
        std::set<const Node*>& nodelist);



    //! Check if region \c id is a region of this simulation environment
    bool contains_region(ID id) const;


    //! Check if an element is an active element of this simulation environment
    /*!
     * \param elem the element to check for
     */
    bool contains_element(const Elem* elem) const;


    //! Check if an ElementSide lies on the boundary of the simulation region
    /*!
     * \deprecated {Use one of is_boundary(), is_outer_boundary() or
     * is_inner_boundary() instead.}
     * This does \em not mean that this side has also a boundary condition
     * assigned.
     */
    bool is_on_boundary(const ElementSide& side) const;


    //! Check if an ElementSide lies on any boundary
    bool is_boundary(const ElementSide& side) const;


    //! Check if an ElementSide lies on a 'real' outer boundary
    bool is_outer_boundary(const ElementSide& side) const;


    //! Check if an ElementSide lies on an inner boundary
    /*!
     * An inner boundary is not really a boundary, but should be
     * considered as an n-1 dimensional domain.
     */
    bool is_inner_boundary(const ElementSide& side) const;


    //! Get the iterator for the first named boundary side
    const BoundarySideIterator boundary_sides_begin(const std::string& name = "") const;


    //! Get the end iterator for the named boundary sides
    const BoundarySideIterator boundary_sides_end(const std::string& name = "") const;


    //! Get the iterator for the first boundary side
    const BoundaryNodeIterator boundary_nodes_begin(void) const;


    //! Get the end iterator for the boundary sides
    const BoundaryNodeIterator boundary_nodes_end(void) const;


    //! Get the iterator for the first boundary
    const BoundaryIterator boundaries_begin(void) const;


    //! Get the end iterator for the boundaries
    const BoundaryIterator boundaries_end(void) const;


    //! Get the iterator for the first boundary element
    /*!
     * \param bd a certain boundary if you want only elements
     * touching this boundary
     */
    const BoundaryElementMap::iterator
      boundary_elements_begin(const Boundary* bd = NULL) const;


    //! Get the end iterator for the boundary elements
    /*!
     * \param bd a certain boundary if you want only elements
     * touching this boundary
     */
    const BoundaryElementMap::iterator
      boundary_elements_end(const Boundary* bd = NULL) const;


    //! Get the iterator for the first boundary element
    /*!
     * \param name a certain boundary name if you want only elements
     * touching this boundary
     */
    const BoundaryElementMap::iterator
      boundary_elements_begin(const std::string& name) const;


    //! Get the end iterator for the boundary elements
    /*!
     * \param name a certain boundary name if you want only elements
     * touching this boundary
     */
    const BoundaryElementMap::iterator
      boundary_elements_end(const std::string& name) const;


    //! Get the iterator for the first element
    const ConstElemIterator elements_begin(void) const;


    //! Get the end iterator for the elements
    const ConstElemIterator elements_end(void) const;


    //! Get the iterator to the first region ID
    const RegionIDIterator region_ids_begin(void) const;


    //! Get the end iterator for the region IDs
    const RegionIDIterator region_ids_end(void) const;


    //! Update the boundary node map
    /*!
     * This method updates the boundary node map from the active elements
     * in the mesh. It should be called after a mesh refinement step.
     */
    void update_boundary_node_map(void);


    //! Create or update the boundary element map
    /*!
     * \param boundaries the boundaries for which the elements should
     * be found
     */
    void update_boundary_element_map(std::set<const Boundary*>
        boundaries = std::set<const Boundary*>());


    //! Update the list of all elements of this simulation
    void update_element_list(void);


    //! Tell if this environment is initialized
    bool is_initialized(void) const;


    //! Tell if this environment is prepared for solve
    bool is_prepared(void) const;


    //! Set the \c unprepared flag
    void invalidate(void);


    //! Get the region numbers
    const std::set<ID>& get_region_ids(void) const;


    //! Get the region names
    void get_region_names(std::set<std::string>& names) const;


    //! Destroy a SimulationEnvironment object
    /*!
     * Checks first if the given pointer is valid yet
     */
    static void destroy(SimulationEnvironment* env);


  private:

    typedef std::set<SimulationEnvironment*> EnvironmentSet;


    //! Disable copy constructor
    SimulationEnvironment(const SimulationEnvironment&);


    //! Disable assignment operator
    SimulationEnvironment& operator=(const SimulationEnvironment&);


    //! Create the list of elements
    void create_element_list(void) TBDLLOCAL;


    //! Creates the boundary element maps
    void create_bc_maps(void) TBDLLOCAL;



    //! Invalidate all relevant environments
    /*!
     * Marks all environments as unprepared
     */
    static void invalidate_all(void) TBDLLOCAL;



    //! The device
    Device* _device;


    //! The mesh
    MeshBase* _mesh;


    //! The region numbers for this simulation
    std::set<ID> _region_numbers;


    //! A set containing all elements belonging to this simulation region
    ElementList _element_list;


    //! A map that assigns boundary ID to boundary condition
    BCMap _bc_map;


    //! The ElementSide for the different boundaries
    /*!
     * Connects ElementSide objects to Boundary pointers
     */
    ElemSideMap _element_side_map;


    //! A map that contains all nodes of the boundaries
    /*!
     * Connects Node pointers to Boundary pointers
     */
    BoundaryNodeMap _node_map;


    //! A map that contains elements touching each boundary
    /*!
     * This map will be created only when calling
     * update_boundary_element_map()
     */
    BoundaryElementMap _bd_elem_map;


    //! Tells if this environment is already initialized
    bool _is_initialized;


    //! Tells if this environment is prepared for solve
    bool _is_prepared;


    //! All environments
    static EnvironmentSet _environments;

};



//
// inline members
//


inline
Device&
SimulationEnvironment::get_device(void)
{
  assert(_device != NULL);
  return *_device;
}

inline
const Device&
SimulationEnvironment::get_device(void) const
{
  assert(_device != NULL);
  return *_device;
}



inline
MeshBase&
SimulationEnvironment::get_mesh(void)
{
  return *_mesh;
}


inline
bool
SimulationEnvironment::is_initialized(void) const
{
  return _is_initialized;
}


inline
bool
SimulationEnvironment::is_prepared(void) const
{
  return _is_prepared;
}


inline
void
SimulationEnvironment::invalidate(void)
{
  _is_prepared = false;
}



inline
bool
SimulationEnvironment::is_on_boundary(const ElementSide& side) const
{
  bool result = true;

  const Elem* neighbour = (side.elem())->neighbor(side.side());

  if (neighbour != NULL)
  {
    ID neighbour_id = neighbour->subdomain_id();
    if (_region_numbers.find(neighbour_id) != _region_numbers.end())
      result = false;
  }

  return result;
}


inline
bool
SimulationEnvironment::is_boundary(const ElementSide& side) const
{
  bool result = false;

  const Elem* neighbour = (side.elem())->neighbor(side.side());

  if ((neighbour == NULL) ||
      (neighbour->subdomain_id() != (side.elem()->subdomain_id())))
    result = true;

  return result;
}



inline
bool
SimulationEnvironment::is_outer_boundary(const ElementSide& side) const
{
  bool result = true;

  const Elem* neighbour = (side.elem())->neighbor(side.side());

  if (neighbour != NULL)
  {
    ID neighbour_id = neighbour->subdomain_id();
    if (_region_numbers.find(neighbour_id) != _region_numbers.end())
      result = false;
  }

  return result;
}


inline
bool
SimulationEnvironment::is_inner_boundary(const ElementSide& side) const
{
  bool result = false;

  const Elem* neighbour = (side.elem())->neighbor(side.side());

  if ((neighbour != NULL) &&
      (neighbour->subdomain_id() != (side.elem()->subdomain_id())))
  {
    ID neighbour_id = neighbour->subdomain_id();
    if (_region_numbers.find(neighbour_id) != _region_numbers.end())
      result = true;
  }

  return result;
}




inline
Boundary*
SimulationEnvironment::get_boundary(const ElementSide& side) const
{
  ElemSideMap::const_iterator it(_element_side_map.find(side));

  if (it != _element_side_map.end())
    return get_boundary(it->second);
  else
    return NULL;
}


inline
Boundary*
SimulationEnvironment::get_boundary(const Node* node) const
{
  std::set<ID> ids;
  _node_map.find_node(node, ids);

  if (!ids.empty())
    return get_boundary(*(ids.begin()));
  else
    return NULL;
}




inline
bool
SimulationEnvironment::contains_region(ID id) const
{
  bool result = true;

  if (_region_numbers.find(id) == _region_numbers.end())
    result = false;

  return result;
}


inline
bool
SimulationEnvironment::contains_element(const Elem* elem) const
{
  bool result = true;

  if (_element_list.find(elem) == _element_list.end())
    result = false;

  return result;
}



/*
inline
const SimulationEnvironment::BoundaryNodeIterator
SimulationEnvironment::boundary_nodes_begin(void) const
{
  return _node_map.begin();
}


inline
const SimulationEnvironment::BoundaryNodeIterator
SimulationEnvironment::boundary_nodes_end(void) const
{
  return _node_map.end();
}
*/


inline
const BoundaryElementMap::iterator
SimulationEnvironment::boundary_elements_begin(const Boundary* bd) const
{
  return _bd_elem_map.elements_begin(bd);
}



inline
const BoundaryElementMap::iterator
SimulationEnvironment::boundary_elements_end(const Boundary* bd) const
{
  return _bd_elem_map.elements_end(bd);
}


inline
const BoundaryElementMap::iterator
SimulationEnvironment::boundary_elements_begin(const std::string& name) const
{
  return boundary_elements_begin(get_boundary(name));
}



inline
const BoundaryElementMap::iterator
SimulationEnvironment::boundary_elements_end(const std::string& name) const
{
  return boundary_elements_end(get_boundary(name));
}


inline
const SimulationEnvironment::RegionIDIterator
SimulationEnvironment::region_ids_begin(void) const
{
  return _region_numbers.begin();
}


inline
const SimulationEnvironment::RegionIDIterator
SimulationEnvironment::region_ids_end(void) const
{
  return _region_numbers.end();
}




inline
const SimulationEnvironment::BoundaryIterator
SimulationEnvironment::boundaries_begin(void) const
{
  return BoundaryIterator(_bc_map.begin());
}


inline
const SimulationEnvironment::BoundaryIterator
SimulationEnvironment::boundaries_end(void) const
{
  return BoundaryIterator(_bc_map.end());
}


inline
const SimulationEnvironment::ConstElemIterator
SimulationEnvironment::elements_begin(void) const
{
  return _element_list.begin();
}


inline
const SimulationEnvironment::ConstElemIterator
SimulationEnvironment::elements_end(void) const
{
  return _element_list.end();
}



inline
void
SimulationEnvironment::update_element_list(void)
{
  create_element_list();
}


inline
void
SimulationEnvironment::get_boundary_nodes(const std::string& name,
    std::set<const Node*>& nodelist)
{
  const Boundary* bd = get_boundary(name);
  get_boundary_nodes(bd, nodelist);
}


inline
const std::set<ID>&
SimulationEnvironment::get_region_ids(void) const
{
  return _region_numbers;
}


// Implementation of BoundarySideIterator
class SimulationEnvironment::BoundarySideIterator
{

  public:
    BoundarySideIterator(void) {}

    BoundarySideIterator(const BoundarySideIterator& it) :
      _iter(it._iter),
      _bdids(it._bdids),
      _mapend(it._mapend) {}

    BoundarySideIterator(const ElemSideMap& emap,
        const ElemSideMap::const_iterator& it,
        const std::set<ID>& ids = std::set<ID>()) :
          _iter(it),
          _bdids(ids),
          _mapend(emap.end())
    {
      while ((_iter != _mapend) && !_bdids.empty() &&
          !_bdids.count(_iter->second))
      {
        ++_iter;
      }
    }

    BoundarySideIterator& operator++(void)
    {
      if (_iter != _mapend)
      {
        ++_iter;
        while ((_iter != _mapend) && !_bdids.empty() &&
            !_bdids.count(_iter->second))
        {
          ++_iter;
        }
      }
      return *this;
    }

    BoundarySideIterator& operator=(const BoundarySideIterator& rhs)
    {
      _iter = rhs._iter;
      _bdids = rhs._bdids;
      _mapend = rhs._mapend;
      return *this;
    }

    bool operator==(const BoundarySideIterator& rhs)
    {
      return ((_iter == rhs._iter) && (_bdids == rhs._bdids));
    }

    bool operator!=(const BoundarySideIterator& rhs)
    {
      return !(*this == rhs);
    }

    const ElementSide& operator*(void)
    {
      return _iter->first;
    }

    const ElemSideMap::const_iterator& operator->(void)
    {
      return _iter;
    }

  private:

    ElemSideMap::const_iterator _iter;

    std::set<ID> _bdids;

    ElemSideMap::const_iterator _mapend;
};

#endif // _SIMULATIONENVIRONMENT_H_
