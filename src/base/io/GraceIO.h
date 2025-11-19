// $Id$


#ifndef _GRACEIO_H_
#define _GRACEIO_H_

#include "tibercad/io/DataOutput.h"
#include "tibercad/base/tiber_dll.h"



//! Write nodal and elemental data using a grace-compatible format
class TBDLLOCAL GraceIO : public DataOutput
{

 public:

  //! Constructor
  GraceIO(void) : DataOutput() {}

  //! Constructor
  /*!
   * \param mesh a reference to a constant mesh object.
   */
  GraceIO(const MeshBase& mesh);


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

   typedef std::map<ID, std::vector<const libMesh::Elem*> > PieceMap;

  void create_pieces(PieceMap& pieces);

};


#endif
