// $Id$


#ifndef _TIBERVTKIO_H_
#define _TIBERVTKIO_H_

#include "TypeDefs.h"
#include "DataOutput.h"
#include "tiber_dll.h"
#include "elem.h"



//! Write data using VTK XML unstructured grid format (.vtu)
class TBDLLOCAL TiberVTKIO : public DataOutput
{

 public:


  //! Constructor
  TiberVTKIO(void) : DataOutput() {}

  //! Constructor
  /*!
   * \param mesh a reference to a constant mesh object.
   */
  TiberVTKIO(const MeshBase& mesh);

  //! Destructor
  ~TiberVTKIO(void) {}



  //! Write a mesh with nodal data
  void write_nodal_data(const std::string& fname,
      const std::vector<double>& soln,
      const std::vector<std::string>& names);


  //! Write a mesh with elemental data
  void write_elemental_data(const std::string& fname,
      const std::vector<double>& soln,
      const std::vector<std::string>& names);



 protected:

  //! The implementation of the writing routine
  virtual void do_write(bool force);



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

  struct VTKElem
  {
    VTKCellType type;
    std::vector<int> connectivity;
  };

  //! Get the VTK cell type for an element
  VTKCellType get_VTK_cell_type(const libMesh::Elem* elem);

  void create_pieces(std::map<ID, std::vector<unsigned int> >& points,
    std::map<ID, std::vector<VTKElem> >& elems);

  //! Write a data array to the given stream
  template <typename T>
  void write_data_array(const std::string& name, int comp, const std::vector<T>& data,
      std::ostream& os);

};


#endif
