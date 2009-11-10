// $Id$


#ifndef _TIBERVTKIO_H_
#define _TIBERVTKIO_H_

#include "mesh_output.h"


// forward declaration
class MeshBase;
class Elem;


//! Write nodal and elemental data using a grace-compatible format 
class TiberVTKIO : public MeshOutput<MeshBase>
{
 public:

    
  //! Constructor
  /*!
   * \param mesh a reference to a constant mesh object.
   */
  TiberVTKIO(const MeshBase& mesh);


  //! Write the mesh to the specified file.
  /*!
   * Does nothing for this format
   */
  virtual void write(const std::string&) {};

  
  //! Write a mesh with nodal data
  virtual void write_nodal_data(const std::string& fname,
      const std::vector<Number>& soln,
      const std::vector<std::string>& names);

  
  //! Write a mesh with elemental data
  virtual void write_elemental_data(const std::string& fname,
      const std::vector<Number>& soln,
      const std::vector<std::string>& names);



 private:

  //! The VTK cell types
  enum VTKCellType
  {
    VTK_UNKNOWN = 0,
    VTK_VERTEX,
    VTK_POLY_VERTEX,
    VTK_LINE,
    VTK_POLY_LINE,
    VTK_TRIANGLE,
    VTK_TRIANGLE_STRIP,
    VTK_POLYGON,
    VTK_PIXEL,
    VTK_QUAD,
    VTK_TETRA,
    VTK_VOXEL,
    VTK_HEXAHEDRON,
    VTK_WEDGE,
    VTK_PYRAMID
  };

  
  //! Get the VTK cell type for an element
  VTKCellType get_VTK_cell_type(const Elem* elem);
  

};

    
#endif
