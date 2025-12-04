/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file BoundaryRegions.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _BOUNDARYREGIONS_H_
#define _BOUNDARYREGIONS_H_

#include "tibercad/geom/MeshRegionInfo.h"
#include "tibercad/geom/ElementSide.h"
#include "tibercad/base/HashMap.h"
#include "tibercad/base/HashSet.h"
#include "tibercad/base/tiber_dll.h"

#include <string>
#include <map>
#include <set>
#include <vector>

namespace libMesh
{
class Elem;
class Node;
}


//! A class containing information on lower dimensional regions
/*!
 * Boundary IDs are assumed to be unique, even if they refer to boundaries
 * of different spatial dimension.
 */
class TBDLLOCAL BoundaryRegions : public MeshRegionInfo
{

  private:
    //! The element side map
    typedef HashMap<ElementSide, ID, ElementSide::hash>::Type ElemSideMap;


  public:

    //! A set of nodes
    typedef HashSet<const libMesh::Node*>::Type NodeSet;

    //! The type of boundary (side, edge node)
    enum Type
    {
        SIDE,    /*!< an element side */
        EDGE,    /*!< an element edge */
        NODE,    /*!< a node */
        INVALID  /*!< invalid */
    };


    //BoundaryRegions(const libMesh::Parallel::Communicator &comm_in) :
    //  MeshRegionInfo(comm_in)
    BoundaryRegions(libMesh::MeshBase &mesh) :
      MeshRegionInfo(mesh)
    {}

    //! Iterator to iterate over a subset of boundaries
    class side_iterator
    {
      public:
        //! Constructor
        side_iterator(ElemSideMap& map, const std::set<ID>& ids,
            bool begin = false) :
          _map(map),
          _ids(ids),
          _iter(_map.begin())
        {
          if (begin && _iter != _map.end())
            while ((_iter != _map.end() && !_ids.empty()
                && !_ids.count(_iter->second)))
              ++_iter;
          else
            _iter = _map.end();
        }

        //! Copy constructor
        side_iterator(const side_iterator& it) :
          _map(it._map),
          _ids(it._ids),
          _iter(it._iter) { }

        //! Prefix increment
        side_iterator& operator++(void)
        {
          do {
            ++_iter;
          }
          while ((_iter != _map.end()) && !_ids.empty()
              && !_ids.count(_iter->second));
          return *this;
        }

        //! Assignment
        side_iterator& operator=(const side_iterator& rhs)
        {
          _map = rhs._map;
          _ids = rhs._ids;
          _iter = rhs._iter;
          return *this;
        }

        //! Comparison
        bool operator==(const side_iterator& rhs)
                       {
          return (_iter == rhs._iter);
                       }

        //! Comparison
        bool operator!=(const side_iterator& rhs)
                       {
          return (_iter != rhs._iter);
                       }

        //! Dereference the iterator
        const ElementSide& operator*(void)
        {
          return _iter->first;
        }

      private:
        ElemSideMap& _map;
        std::set<ID> _ids;
        ElemSideMap::iterator _iter;
    };


    //! Destructor
    ~BoundaryRegions(void);


    //! Add an element side with boundary ID \c id
    void add_side(const libMesh::Elem* elem, unsigned int side, ID id);


    //! Add an element edge with boundary ID \c id
    void add_edge(const libMesh::Elem* elem, unsigned int edge, ID id);


    //! Add a node with boundary ID \c id
    void add_node(const libMesh::Node* node, ID id);


    //! Get the ID associated to a side
    ID get_side_id(const libMesh::Elem* elem, unsigned int side) const;


    //! Get contiguous regions for a given side ID
    /*!
     * An empty set indicates that boundary \c id does not exist.
     * If one of the IDs is INVALID_ID it means that \c ID belongs
     * to an external boundary.
     */
    const std::set<ID>& get_contiguous_regions_for_side(ID id);


    //! Get the ID associated to a edge
    ID get_edge_id(const libMesh::Elem* elem, unsigned int edge) const;


    //! Get the ID associated to a node
    ID get_node_id(const libMesh::Node* node) const;


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


    void get_bc_node_map(std::map<ID, std::vector<unsigned int> >& nodemap) const;


    //! Get the side iterator for a given set of boundary IDs
    side_iterator sides_begin(const std::set<ID>& ids);

    //! Get the past-the-end iterator for a given set of boundary IDs
    side_iterator sides_end(const std::set<ID>& ids);


  protected:

    virtual void do_broadcast(void);


  private:

    //! The subdomain IDs on both sides of a boundary
    /*!
     * We use std::set as value type as this can be easily compared
     */
    typedef HashMap<ID, std::set<ID> >::Type BDMatMap;

    //! The element edge map
    typedef HashMap<ElementEdge, ID, ElementEdge::hash>::Type ElemEdgeMap;

    //! The node map
    typedef HashMap<const libMesh::Node*, ID>::Type NodeMap;


    //! The sides
    ElemSideMap _sides;


    //! All IDs associated to sides
    std::set<ID> _side_ids;


    //! The regions touching each boundary region
    BDMatMap _contiguous_regions;


    //! The edges
    ElemEdgeMap _edges;


    //! All IDs associated to edges
    std::set<ID> _edge_ids;


    //! The node map
    NodeMap _nodes;


    //! All IDs associated to nodes
    std::set<ID> _node_ids;


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
BoundaryRegions::get_side_id(const libMesh::Elem* elem, unsigned int side) const
{
  ID id = INVALID_ID;
  ElemSideMap::const_iterator it(_sides.find(ElementSide(elem, side)));
  if (it != _sides.end())
    id = it->second;

  return id;
}


inline
ID
BoundaryRegions::get_edge_id(const libMesh::Elem* elem, unsigned int edge) const
{
  ID id = INVALID_ID;
  ElemEdgeMap::const_iterator it(_edges.find(ElementEdge(elem, edge)));
  if (it != _edges.end())
    id = it->second;

  return id;
}


inline
ID
BoundaryRegions::get_node_id(const libMesh::Node* node) const
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


inline
BoundaryRegions::side_iterator
BoundaryRegions::sides_begin(const std::set<ID>& ids = std::set<ID>())
{
  return side_iterator(_sides, ids, true);
}

inline
BoundaryRegions::side_iterator
BoundaryRegions::sides_end(const std::set<ID>& ids = std::set<ID>())
{
  return side_iterator(_sides, ids);
}


#endif // _BOUNDARYREGIONS_H_
