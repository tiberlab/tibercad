// $Id$

#ifndef _BOUNDARYREGIONS_H_
#define _BOUNDARYREGIONS_H_

#include "TypeDefs.h"
#include "ElementSide.h"
#include "HashMap.h"
#include "HashSet.h"

#include <string>
#include <map>
#include <vector>

class Elem;
class Node;


//! A class containing information on lower dimensional regions
/*!
 * Boundary IDs are assumed to be unique, even if they refer to boundaries
 * of different spatial dimension.
 */
class BoundaryRegions
{

  public:

    //! A set of nodes
    typedef TiberCad::HashSet<const Node*>::Type NodeSet;


    //! Destructor
    ~BoundaryRegions(void);


    //! Add a element side with boundary ID \c id
    void add_side(const Elem* elem, unsigned int side, ID id);


    //! Add a element edge with boundary ID \c id
    void add_edge(const Elem* elem, unsigned int edge, ID id);


    //! Add a node with boundary ID \c id
    void add_node(const Node* node, ID id);


    //! Add a name to a boundary ID
    /*!
     * Nonexistent IDs are not added.
     */
    void set_name(ID id, const std::string& name);


    //! Get the set of all side nodes of a certain boundary ID
    void get_side_nodes(ID id, NodeSet& nodes) const;


    //! Get the set of all edge nodes of a certain boundary ID
    void get_edge_nodes(ID id, NodeSet& nodes) const;


    //! Get the set of all nodes of a certain node-type boundary ID
    void get_nodes(ID id, NodeSet& nodes) const;


    //! Get the ID for a boundary name
    ID get_id(const std::string& name) const;


    //! Cleanup the structure
    void clear(void);


    void get_bc_node_map(std::map<unsigned int, std::vector<ID> >& nodemap) const;


  private:

    //! The element side map
    typedef TiberCad::HashMap<ElementSide, ID, ElementSide::hash>::Type ElemSideMap;

    //! The element edge map
    typedef TiberCad::HashMap<ElementEdge, ID, ElementEdge::hash>::Type ElemEdgeMap;

    //! The node map
    typedef TiberCad::HashMap<const Node*, ID>::Type NodeMap;

    //! ID to name map
    typedef TiberCad::HashMap<ID, std::string>::Type IDToNameMap;

    //! name to ID map
    typedef TiberCad::HashMap<std::string, ID>::Type NameToIDMap;


    //! The sides
    ElemSideMap _sides;


    //! The edges
    ElemEdgeMap _edges;


    //! The node map
    NodeMap _nodes;


    //! The names corresponding to the boundary IDs
    /*!
     * This map contains all IDs, even those without a name
     */
    IDToNameMap _ids_to_names;

    //! The IDs corresponding to the names
    NameToIDMap _names_to_ids;

};


//
// inline members
//



inline
BoundaryRegions::~BoundaryRegions(void)
{
  clear();
}

#endif // _BOUNDARYREGIONS_H_
