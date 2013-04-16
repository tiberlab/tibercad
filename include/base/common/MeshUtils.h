// $Id$


#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

#include "TypeDefs.h"
#include "TensorGrid.h"
#include "point.h"
#include "tiber_dll.h"

#include "auto_ptr.h"

#include <set>
#include <vector>
#include <map>

class Elem;
class MeshBase;


//! A few utilities for operations on a mesh
class MeshUtils
{

  public:


    //! A helper class to map grids onto a TensorGrid
    class GridMapper
    {
      public:

        //! Obtain the GridMapper object for a given mesh
        static GridMapper& get_mapper(const MeshBase* mesh);

        //! Get the element a given point is in
        const Elem* get_element(const Point& point) const;


      private:

        typedef std::vector<std::vector<const Elem*>>  ElementList;

        //! Constructor
        GridMapper(const MeshBase* mesh);

        //! The real mesh
        const MeshBase* _mesh;

        //! The tensor grid
        TensorGrid _tensor_grid;

        //! The datastructure containing the element lists
        ElementList _elem_list;

        //! Setup the mapping
        void setup(void);


        //! A static list of all GridMapper objects
        static std::map<const MeshBase*, GridMapper> _mappers;

    };


    //! Get all the subdomain IDs present in the mesh
    /*!
     * This method iterates over all mesh elements and puts the subdomain id
     * into the provided set.
     */
    static void get_subdomain_ids(MeshBase& mesh, std::set<ID>& subdomain_ids);


    //! Preliminary check to see if a point can belong to an element
    /*!
     * If true, the point can belong to the element (but not for sure).
     * If false, the point does not belong to element.
     * It's much faster then exact check. As it uses a parallepipedal box
     * around the element, it only performs a simple matrix-vector product
     */
    static bool may_belong_to_element(const Elem* element, const Point& point);


    //! Get the element a given point lies in
    static const Elem* search_element(const MeshBase* mesh, const Point& point);


    //! Get the outer normal on an element side
    static Point get_outer_normal(const Elem* elem, int side);

    //! Create the boundary mesh for a given volume mesh
    static AutoPtr<MeshBase> create_boundary_mesh(const MeshBase& mesh);


  private:


    //! This class contains only static methods
    MeshUtils(void);


};


#endif // _MESHUTILS_H_
