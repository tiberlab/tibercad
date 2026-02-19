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
 * \file ReadComsol.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */


#ifndef TC_READCOMSOL_H
#define TC_READCOMSOL_H

#include "tibercad/base/tiber_dll.h"

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
class TC_DLLOCAL ReadComsol : public libMesh::MeshInput<libMesh::MeshBase>

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


#endif // TC_READCOMSOL_H
