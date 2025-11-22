// $Id$

#ifndef __tecplot_io_cell_h__
#define __tecplot_io_cell_h__

// Local includes
//#include "libmesh_common.h"
//#include "mesh_output.h"
#include "tecplot_io.h"
#include "tibercad/base/tiber_dll.h"

#include "tibercad/base/libMeshDefs.h"




/**
 * This class implements writing meshes with cell data in the Tecplot format.
 * (derived from TecplotIO)
 *
 * @author Benjamin S. Kirk, 2004, Optolab  2006
 * Optolab  2006
 */

// ------------------------------------------------------------
// TecplotIO_cell class definition
//class TecplotIO_cell : public MeshOutput<MeshBase>
class TBDLLOCAL TecplotIO_cell : public libMesh::TecplotIO

{
 public:

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   * The optional parameter \p binary can be used to switch
   * between ASCII (\p false, the default) or binary (\p true)
   * output files.
   */
  // TecplotIO_cell (const MeshBase&, const bool binary=false); 
  TecplotIO_cell (const libMesh::MeshBase& mesh, const bool binary=false)
    : libMesh::TecplotIO(mesh, binary) {} ;


  void write_cell_data (const std::string& fname,
				  const std::vector<Number>& soln,
				 const std::vector<std::string>& names);

  void write_ascii_cell (const std::string& fname,
			     const std::vector<Number>* v,
				  const std::vector<std::string>* solution_names);

};



#endif // #define __tecplot_io_cell_h__
