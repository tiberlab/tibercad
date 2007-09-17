// $Id$


#ifndef _VTKIO_H_
#define _VTKIO_H_

#include "mesh_output.h"


// forward declaration
class MeshBase;


//! Write nodal and elemental data using a grace-compatible format 
class VTKIO : public MeshOutput<MeshBase>
{
 public:

    
  //! Constructor
  /*!
   * \param mesh a reference to a constant mesh object.
   */
  VTKIO(const MeshBase& mesh);


  //! Write the mesh to the specified file.
  /*!
   * Does nothing for this format
   */
  virtual void write(const std::string&) {};

  
  //! Write a mesh with nodal data
  virtual void write_nodal_data(const std::string&,
      const std::vector<Number>&,
      const std::vector<std::string>&);

  
  //! Write a mesh with elemental data
  virtual void write_elemental_data(const std::string&,
      const std::vector<Number>&,
      const std::vector<std::string>&);



 private:

};

    
#endif
