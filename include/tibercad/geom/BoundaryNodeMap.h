// $Id$

#ifndef _BOUNDARYNODEMAP_H_
#define _BOUNDARYNODEMAP_H_

#include "tibercad/base/TypeDefs.h"

#include <map>
#include <set>
#include <vector>

namespace libMesh
{
class Node;
}


//! Contains a map connecting boundary IDs to set of nodes
class BoundaryNodeMap
{

  public:

    //! A typedef for convenience
    typedef std::set<const libMesh::Node*> NodeSet;



    //! Constructor
    BoundaryNodeMap(void);

    //! Get the set of nodes for a certain boundary ID
    const NodeSet& get_nodes(ID id) const;

    //! Add a node set
    //void add_node_set(ID id, const NodeSet& node_set);

    //! Add a single node
    void add_node(ID id, const libMesh::Node* node);

    //! Find the boundaries corresponding to a node
    /*!
     * \return true if the node is on a boundary
     */
    bool find_node(const libMesh::Node* node, std::set<ID>& ids) const;

    //! Find the boundaries corresponding to a node
    /*!
     * \return true if the node is on a boundary
     */
    bool find_node(const libMesh::Node* node, std::vector<ID>& ids) const;




  private:

    std::map<ID, NodeSet> _map;
};


//
// inline methods
//

inline
void
BoundaryNodeMap::add_node(ID id, const libMesh::Node* node)
{
  _map[id].insert(node);
}


#endif // _BOUNDARYNODEMAP_H_
