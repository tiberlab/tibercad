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
 * \file MeshReader.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/io/MeshReader.h"
#include "Read_ISE_grid/ReadISEGrid.h"
#include "ReadGMSH.h"
#include "ReadComsol.h"
#include "tibercad/geom/BoundaryRegions.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/utils/Utils.h"
#include "tibercad/io/Messages.h"

#include "mesh_base.h"
#include "mesh_communication.h"

#include <sstream>

using namespace std;




void MeshReader::read_mesh(const string& filename, libMesh::MeshBase& mesh,
    MeshRegionInfo& region_info, BoundaryRegions& bd_regions)
{

  Utils::Timer tt;

  // we read only on processor 0
  if (mesh.comm().rank() == 0)
  {
    if  (filename.rfind(".grd") < filename.size())
    {
      ReadISEGrid ise_mesh(mesh, region_info, bd_regions);
      ise_mesh.read(filename);
    }
    else if (filename.rfind(".msh") < filename.size())
    {
      ReadGMSH msh_mesh(mesh, region_info, bd_regions);
      msh_mesh.read(filename);
    }
    else if (filename.rfind(".mphtxt") < filename.size())
    {
      ReadComsol comsol_mesh(mesh, region_info, bd_regions);
      comsol_mesh.read(filename);
    }
    else
    {
      ostringstream os;
      os << filename << " has unknown mesh file format." << endl;
      os << "Known formats are: ISE grid (.grd), GMSH (.msh)";
      throw InitFailedException(os.str());
    }

    ostringstream os;
    os << "mesh reading time: " << tt.elapsed_string();
    Messages::info(os.str());

  }

  tt.reset();
  // broadcast it to the other processors
  libMesh::MeshCommunication().broadcast(mesh);
  region_info.broadcast();
  bd_regions.broadcast();

  if (mesh.comm().size() > 1)
  {
    ostringstream os;
    os << "mesh broadcast time: " << tt.elapsed_string();
    Messages::info(os.str());
  }
  tt.reset();

  // now prepare it for use
  mesh.prepare_for_use();
  {
    ostringstream os;
    os << "mesh preparation time: " << tt.elapsed_string();
    Messages::info(os.str());
  }
  tt.reset();

  bd_regions.prepare_for_use();
  {
    ostringstream os;
    os << "bd_regions preparation time: " << tt.elapsed_string();
    Messages::info(os.str());
  }
}
