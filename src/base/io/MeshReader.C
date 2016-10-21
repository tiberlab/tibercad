// $Id$

#include "MeshReader.h"
#include "ReadISEGrid.h"
#include "ReadGMSH.h"
#include "ReadComsol.h"
#include "BoundaryRegions.h"
#include "InitFailedException.h"

#include "mesh_base.h"
#include "mesh_communication.h"

#include <sstream>

using namespace std;




void MeshReader::read_mesh(const string& filename, libMesh::MeshBase& mesh,
    MeshRegionInfo& region_info, BoundaryRegions& bd_regions)
{

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
  }

  // broadcast it to the other processors
  libMesh::MeshCommunication().broadcast(mesh);
  region_info.broadcast();
  bd_regions.broadcast();

  // now prepare it for use
  mesh.prepare_for_use();

  bd_regions.prepare_for_use();
}
