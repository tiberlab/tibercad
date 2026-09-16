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
 * \file MeshUtils.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_MESHUTILS_H
#define TC_MESHUTILS_H

#include "tibercad/base/TypeDefs.h"
#include "tibercad/geom/TensorGrid.h"
#include "tibercad/base/tiber_dll.h"

#include "auto_ptr.h"
#include "elem.h"
#include "mesh_base.h"
#include "point.h"

#include <set>
#include <vector>
#include <map>





//! A few utilities for operations on a mesh
class MeshUtils
{

  public:

    class GridMapper;





    //! Get all the subdomain IDs present in the mesh
    /*!
     * This method iterates over all mesh elements and puts the subdomain id
     * into the provided set.
     */
    static void get_subdomain_ids(libMesh::MeshBase& mesh, std::set<ID>& subdomain_ids);


    //! Preliminary check to see if a point can belong to an element
    /*!
     * If true, the point can belong to the element (but not for sure).
     * If false, the point does not belong to element.
     * It's much faster then exact check. As it uses a parallepipedal box
     * around the element, it only performs a simple matrix-vector product
     */
    static bool may_belong_to_element(const libMesh::Elem* element, libMesh::Point& point);


    //! Get the element a given point lies in
    static const libMesh::Elem* search_element(const libMesh::MeshBase* mesh, const libMesh::Point& point);


    //! Get the outer normal on an element side
    static libMesh::Point get_outer_normal(const libMesh::Elem* elem, int side);

    //! Create the boundary mesh for a given volume mesh
    static std::unique_ptr<libMesh::MeshBase> create_boundary_mesh(const libMesh::MeshBase& mesh);


    //! Triangulate a given point set
    /*!
     * This creates a mesh containing a triangulation of a given point set in 2D,
     * without adding additional points. It follows the description in
     * http://paulbourke.net/papers/triangulate/
     *
     * \param mesh a mesh object containing only the points
     */
    static void triangulate_point_set(libMesh::MeshBase& mesh);

  private:


    //! This class contains only static methods
    MeshUtils(void);



  public:

    //! A helper class to map grids onto a TensorGrid
    /*!
     * It is assumed that there are no overlapping elements.
     * In presence of e.g. automatically generated quantum contacts this may
     * not be the case, but such pieces of the mesh should append elements to
     * the end, which would not be seen normally by the GridMapper.
     * If the GridMapper should include such regions, the set of IDs has to be
     * specified explicitly.
     */
    class GridMapper
    {
      public:

        //! Destructor
        ~GridMapper(void);

        //! Obtain the GridMapper object for a given mesh
        /*!
         * \param regions a subset of the region IDs of \c mesh
         */
        static GridMapper& get_mapper(const libMesh::MeshBase* mesh,
            const std::set<ID>& regions = std::set<ID>());

        //! Obtain the GridMapper object for a given mesh
        /*!
         * \param regions a subset of the region IDs of \c mesh
         */
        static GridMapper& get_mapper(const libMesh::MeshBase& mesh,
            const std::set<ID>& regions = std::set<ID>());

        //! Get the element a given point is in
        const libMesh::Elem* get_element(const libMesh::Point& point) const;


      private:

        typedef std::vector<std::vector<const libMesh::Elem*>>  ElementList;

        //! Default constructor is disabled
        GridMapper(void);

        //! Constructor
        explicit GridMapper(const libMesh::MeshBase* mesh,
            const std::set<ID>& regions = std::set<ID>());

        //! The real mesh
        const libMesh::MeshBase* _mesh;

        //! The subset of mesh region IDs
        std::set<ID> _regids;

        //! The tensor grid
        TensorGrid _tensor_grid;

        //! The datastructure containing the element lists
        ElementList _elem_list;

        //! Setup the mapping
        void setup(void);


        //! A static list of all GridMapper objects
        static std::multimap<const libMesh::MeshBase*, GridMapper*> _mappers;

    };
};


#endif // TC_MESHUTILS_H
