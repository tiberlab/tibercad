// $Id$


#ifndef _MESHUTILS_H_
#define _MESHUTILS_H_

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

  private:

    //! This class contains only static methods
    MeshUtils(void);

};


#endif // _MESHUTILS_H_
