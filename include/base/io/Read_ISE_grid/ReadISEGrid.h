#ifndef _READISEGRID_H_
#define _READISEGRID_H_


#include <iomanip>
#include <set>
#include <vector>
#include <map>

#include <iostream> //for I/O interaction
#include <fstream>  //for file streaming
#include <InitFailedException.h>

#include "TypeDefs.h"
#include "tiber_dll.h"

//  LibMesh  include files
#include "mesh_input.h"

class MeshBase;
class MeshRegionInfo;
class BoundaryRegions;
class ISE_Element;
class ISE_Vertex;
class ISE_Edge;
class ISE_Face;

//! Reads an ISE mesh (*.grd file).
class TBDLLOCAL ReadISEGrid : public MeshInput<MeshBase>
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
  MeshInput<MeshBase>(mesh),
  _reg_info(reg_info),
  _bd_regions(bd_regions)
{}




#endif /* _READISEGRID_H_ */
