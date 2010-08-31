// $Id$


#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

#include "TypeDefs.h"
#include "point.h"
#include "tiber_dll.h"

#include <set>

class Elem;
class MeshBase;


//! A few utilities for operations on a mesh
class TBDLLOCAL MeshUtils
{

  public:

    //! Get all the subdomain IDs present in the mesh
    /*!
     * This method iterates over all mesh elements and puts the subdomain id
     * into the provided set.
     */
    static void get_subdomain_ids(MeshBase& mesh, std::set<ID>& subdomain_ids);

    //!Preliminary check to see if a point can belong to an element
    /*!
     * If true, the point can belong to the element (but not for sure).
     * If false, the point does not belong to element.
     * It's much faster then exact check. As it uses a parallepipedal box
     * around the element, it only performs a simple matrix-vector product
     */
    static bool may_belong_to_element(const Elem* element, Point& point);

  private:

    //! This class contains only static methods
    MeshUtils(void);

};


#endif // _MESHUTILS_H_
