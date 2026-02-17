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
 * \file ReadISEGrid.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef TC_READISEGRID_H
#define TC_READISEGRID_H

#include "tibercad/base/InitFailedException.h"
#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/libMeshDefs.h"
#include "tibercad/base/tiber_dll.h"

//  LibMesh  include files
#include "libmesh/mesh_input.h"

#include <iomanip>
#include <set>
#include <vector>
#include <map>

#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming

class MeshRegionInfo;
class BoundaryRegions;
class ISE_Element;
class ISE_Vertex;
class ISE_Edge;
class ISE_Face;



//! Reads an ISE mesh (*.grd file).
class TBDLLOCAL ReadISEGrid : public libMesh::MeshInput<MeshBase>
{

 public:

  /*!
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  ReadISEGrid(MeshBase& mesh, MeshRegionInfo& reg_info,
      BoundaryRegions& bd_regions);


  //! Virtual Destructor.
  /*!
    Deallocate pointers.
  */
  virtual ~ReadISEGrid() {};

  /*!
   * Reads in a mesh in the Gmsh *.msh format
   * from the ASCII file given by name.
   */
  virtual void read(const std::string& name);



 private:

  /*!
   *  Integrity Check. Controls if file Version and Type are correct.
   *  Otherwise aborts program.
   */
  void integrity_check(float ver, std::string tp);

  /*!
    Reads Input file.
  */
  void scan_grid_file(std::istream& ISE_INPUT);



  //! Mesh region info
  MeshRegionInfo& _reg_info;


  //! The object to hold boundary region information
  BoundaryRegions& _bd_regions;


};


inline
ReadISEGrid::ReadISEGrid(MeshBase& mesh, MeshRegionInfo& reg_info,
    BoundaryRegions& bd_regions) :
  libMesh::MeshInput<MeshBase>(mesh),
  _reg_info(reg_info),
  _bd_regions(bd_regions)
{}




#endif /* _READISEGRID_H_ */
