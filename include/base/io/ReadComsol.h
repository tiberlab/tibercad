// $Id: ReadComsol.h 2057 2010-08-31 13:00:49Z maufder $

#ifndef __READCOMSOL_H_
#define __READCOMSOL_H_

#include "tiber_dll.h"

// libMesh includes
#include "libmesh_common.h"
#include "mesh_input.h"

// Forward declarations
class BoundaryRegions;
class MeshRegionInfo;

namespace libMesh
{
  class MeshBase;
}



/*!
 * This class implements reading of meshes in the Comsol format.
 */
class TBDLLOCAL ReadComsol : public libMesh::MeshInput<libMesh::MeshBase>

{
  public:

    /*!
     * Constructor.  Takes a non-const Mesh reference which it
     * will fill up with elements via the read() command.
     */
    ReadComsol(libMesh::MeshBase& mesh, MeshRegionInfo& reg_info, BoundaryRegions& bd_regions);

    /*!
     * Reads in a mesh in the Comsol *.mphtxt format
     * from the ASCII file given by name.
     */
    virtual void read(const std::string& name);


  private:

    /*!
     * Implementation of the read() function.  This function
     * is called by the public interface function and implements
     * reading the file.
     */
    void read_mesh(std::istream& in);


    //! Mesh region info
    MeshRegionInfo& _reg_info;


    //! The object to hold boundary region information
    BoundaryRegions& _bd_regions;


};


inline
ReadComsol::ReadComsol(libMesh::MeshBase& mesh, MeshRegionInfo& reg_info,
    BoundaryRegions& bd_regions) :
  MeshInput<libMesh::MeshBase>(mesh),
  _reg_info(reg_info),
  _bd_regions(bd_regions)
{}


#endif // __READCOMSOL_H_
