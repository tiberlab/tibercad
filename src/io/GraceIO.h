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
 * \file GraceIO.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */



#ifndef TC_GRACEIO_H
#define TC_GRACEIO_H

#include "tibercad/io/DataOutput.h"
#include "tibercad/base/tiber_dll.h"



//! Write nodal and elemental data using a grace-compatible format
class TC_DLLOCAL GraceIO : public DataOutput
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
