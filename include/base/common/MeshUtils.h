// $Id$


#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

#include "TypeDefs.h"

#include <set>

class Mesh;
class MeshData_elements;


//! A few utilities for operations on a mesh
class MeshUtils
{

  public:

    //! Assign the right subdomain id to every element
    /*!
     * This method loops over all elements in \c meshdata and assigns
     * the subdomain id to it.
     *
     * \param mesh the mesh itself
     * \param meshdata the meshdata object
     */
    static void assign_subdomain_ids(Mesh& mesh, MeshData_elements& meshdata);


    //! Get all the subdomain IDs present in the mesh
    /*!
     * This method iterates over all mesh elements and puts the subdomain id
     * into the provided set.
     */
    static void get_subdomain_ids(Mesh& mesh, std::set<ID>& subdomain_ids);

  private:

    //! This class contains only static methods
    MeshUtils(void);

};


#endif // _MESHUTILS_H_
