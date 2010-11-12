// $Id$

#ifndef _BOUNDARYREGIONS_H_
#define _BOUNDARYREGIONS_H_

#include "MeshRegionInfo.h"
#include "ElementSide.h"
#include "HashMap.h"
#include "HashSet.h"
#include "tiber_dll.h"

#include <string>
#include <map>
#include <set>
#include <vector>

class Elem;
class Node;


//! A class containing information on lower dimensional regions
/*!
 * Boundary IDs are assumed to be unique, even if they refer to boundaries
 * of different spatial dimension.
 */
class TBDLLOCAL BoundaryRegions : public MeshRegionInfo
{

  public:

    //! A set of nodes
    typedef HashSet<const Node*>::Type NodeSet;

    //! The type of boundary (side, edge node)
    enum Type
    {
        SIDE,    /*!< an element side */
        EDGE,    /*!< an element edge */
        NODE,    /*!< a node */
        INVALID  /*!< invalid */
    };


    //! Destructor
    ~BoundaryRegions(void);


    //! Add an element side with boundary ID \c id
    void add_side(const Elem* elem, unsigned int side, ID id);


    //! Add an element edge with boundary ID \c id
    void add_edge(const Elem* elem, unsigned int edge, ID id);


    //! Add a node with boundary ID \c id
    void add_node(const Node* node, ID id);


    //! Get the ID associated to a side
    ID get_side_id(const Elem* elem, unsigned int side) const;


    //! Get contiguous regions for a given side ID
    /*!
     * The size of the set returned is always 2 (unless the side ID does
     * not exist, in which case an empty set is returned).
     * If one of the IDs is * INVALID_ID it means that \c ID belongs
     * to an external boundary.
     */
    const std::set<ID>& get_contiguous_regions_for_side(ID id);


    //! Get the ID associated to a edge
    ID get_edge_id(const Elem* elem, unsigned int edge) const;


    //! Get the ID associated to a node
    ID get_node_id(const Node* node) const;


    //! Get the set of all side nodes of a certain boundary ID
    void get_side_nodes(ID id, NodeSet& nodes) const;


    //! Get the set of all edge nodes of a certain boundary ID
    void get_edge_nodes(ID id, NodeSet& nodes) const;


    //! Get the set of all nodes of a certain node-type boundary ID
    void get_nodes(ID id, NodeSet& nodes) const;


    //! Check if an ID is associated to a side
    bool is_side(ID id) const;


    //! Check if an ID is associated to a edge
    bool is_edge(ID id) const;


    //! Check if an ID is associated to a node
    bool is_node(ID id) const;


    //! Cleanup the structure
    void clear(void);


    //! Prepare the structure for use
    /*!
     * Reassign side IDs such that each boundary between different regions
     * has a unique ID.
     * This method has to be called \em after set_name() and \em before using the
     * object for identifying boundary models.
     */
    void prepare_for_use(void);


    void get_bc_node_map(std::map<unsigned int, std::vector<ID> >& nodemap) const;


  private:

    //! The element side map
    typedef HashMap<ElementSide, ID, ElementSide::hash>::Type ElemSideMap;

    //! The subdomain IDs on both sides of a boundary
    /*!
     * We use std::set as value type as this can be easily compared
     */
    typedef HashMap<ID, std::set<ID> >::Type BDMatMap;

    //! The element edge map
    typedef HashMap<ElementEdge, ID, ElementEdge::hash>::Type ElemEdgeMap;

    //! The node map
    typedef HashMap<const Node*, ID>::Type NodeMap;


    //! The sides
    ElemSideMap _sides;


    //! All IDs associated to sides
    IDHashSet _side_ids;


    //! The regions touching each boundary region
    BDMatMap _contiguous_regions;


    //! The edges
    ElemEdgeMap _edges;


    //! All IDs associated to sides
    IDHashSet _edge_ids;


    //! The node map
    NodeMap _nodes;


    //! All IDs associated to sides
    IDHashSet _node_ids;


};


//
// inline members
//


inline
bool
BoundaryRegions::is_side(ID id) const
{
  return (_side_ids.find(id) != _side_ids.end());
}

inline
bool
BoundaryRegions::is_edge(ID id) const
{
  return (_edge_ids.find(id) != _edge_ids.end());
}

inline
bool
BoundaryRegions::is_node(ID id) const
{
  return (_node_ids.find(id) != _node_ids.end());
}



inline
ID
BoundaryRegions::get_side_id(const Elem* elem, unsigned int side) const
{
  ID id = INVALID_ID;
  ElemSideMap::const_iterator it(_sides.find(ElementSide(elem, side)));
  if (it != _sides.end())
    id = it->second;

  return id;
}


inline
ID
BoundaryRegions::get_edge_id(const Elem* elem, unsigned int edge) const
{
  ID id = INVALID_ID;
  ElemEdgeMap::const_iterator it(_edges.find(ElementEdge(elem, edge)));
  if (it != _edges.end())
    id = it->second;

  return id;
}


inline
ID
BoundaryRegions::get_node_id(const Node* node) const
{
  ID id = INVALID_ID;
  NodeMap::const_iterator it(_nodes.find(node));
  if (it != _nodes.end())
    id = it->second;

  return id;
}


inline
const std::set<ID>&
BoundaryRegions::get_contiguous_regions_for_side(ID id)
{
  return _contiguous_regions[id];
}

inline
BoundaryRegions::~BoundaryRegions(void)
{
  clear();
}

#endif // _BOUNDARYREGIONS_H_
