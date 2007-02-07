// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef __mesh_data_elements_h__
#define __mesh_data_elements_h__

/* #ifndef __mesh_data_h__ */
/* #define __mesh_data_h__ */

// C++ includes
#include <map>
#include <vector>

// Local Includes
#include "libmesh.h"
#include "enum_xdr_mode.h"
#include "mesh_data.h"



// Forward Declarations
/* class Node; */
/* class Elem; */
class MeshBase;
/* class UNVIO; */
/* class TetGenIO; */
/* class MeshDataUnvHeader; */
class MeshData;
/* class BoundaryInfo; */

class MeshData_elements;







/**
 * Yet another Mesh-something class...  What's this good for: 
 * \p MeshData handles actual data and the corresponding I/O
 * on entities (nodes, elements) of meshes.  
 * \p MeshBase owns a \p MeshData for dealing with files
 * that contain nodal or element-oriented data, numbered in the same 
 * format as the corresponding mesh file (when activated) or with
 * the \p MeshBase element and node ids (when in compatibility mode).  
 * To use \p MeshData, it has to be activated or the compatibility 
 * mode has to be enabled.
 *
 * @author Daniel Dreyer, 2003
 */

// ------------------------------------------------------------
// MeshData class definition
//class MeshData
class MeshData_elements : public MeshData { 

public:


  MeshData_elements  (const MeshBase& m): MeshData(m){};  // ok :  constructor derived from MeshData constructor !!


  /**
   * Destructor.
   */
  //  ~MeshData ();
  ~MeshData_elements ();   //   destructor MUST  be  implemented in  derived class (see mesh_data_elements.C)

   



  //---------------------------------------------
  // NEW :
  // translate  for element data


  void translate_elem_data (const MeshBase& out_mesh,
			    std::vector<Number>& data_values,
			    std::vector<std::string>& data_names) const;

  //------------------------------------------------------------------------


 //----------------------------------------------------------
  // Element-associated data
  /**
   * @returns the \f$ i^{th} \f$ value (defaults to 0) associated 
   * with element \p elem.  Returns \p libMesh::zero when there
   * is no data for \p elem in the map.
   */
  Number operator() (const Elem* elem, 
		     const unsigned int i=0) const;









};


//-------------------------------------------------------------
// element data inline methods
inline
Number MeshData_elements::operator() (const Elem* elem, 
			     const unsigned int i) const
{
  assert (_active || _compatibility_mode);
  assert (_elem_data_closed);

  std::map<const Elem*, 
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

/*  std::map<const Elem*,  */
/*            std::vector<Number> >::const_iterator pos = _elem_data.find(elem ->top_parent() ); */

 //  for mesh refinements :  find  elem_data  of  top parent (root) of  current  element !!!
  //  does not  work:  why  ????
 //  ??????????????????????????????????????????????//

  if (pos == _elem_data.end())
    return libMesh::zero;
  
  // we only get here when pos != _elem_data.end()  
  assert (i < pos->second.size());
  return pos->second[i];
}





#endif
