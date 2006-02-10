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


/*   //---------------------------------------------------------- */
/*   // convenient typedefs */
/*   /** */
/*    * A const iterator over the nodal data entries of */
/*    * \p MeshData.  Use this when a loop over all \p Node* */
/*    * in the \p MeshData is wanted.  Note that only const versions  */
/*    * are provided.  Also these iterators should @e not be  */
/*    * confused with the \p node_iterators provided */
/*    * for the \p Mesh classes! */
/*    */ 
/*   typedef std::map<const Node*, std::vector<Number> >::const_iterator const_node_data_iterator; */

/*   /** */
/*    * A const iterator over the element-associated data entries of */
/*    * \p MeshData.  Use this when a loop over all \p Node* */
/*    * in the \p MeshData is wanted.  Note that only const versions  */
/*    * are provided.  Also these iterators should @e not be  */
/*    * confused with the \p node_iterators provided */
/*    * for the \p Mesh classes! */
/*    */ 
/*   typedef std::map<const Elem*, std::vector<Number> >::const_iterator const_elem_data_iterator; */


  //----------------------------------------------------------
  /**
   * Default Constructor.  Takes const reference
   * to the mesh it belongs to.
   */
  //  MeshData (const MeshBase& m);

 

  MeshData_elements  (const MeshBase& m): MeshData(m){};  // ok :  constructor derived from MeshData constructor !!


  /**
   * Destructor.
   */
  //  ~MeshData ();
  ~MeshData_elements ();   //   destructor MUST  be  implemented in  derived class (see mesh_data_elements.C)

   



/*   /** */
/*    * When \p MeshData should be used, it has to be activated */
/*    * first, @e prior to reading in a mesh with the \p Mesh::read() */
/*    * methods.  Optionally takes a string that should help the user */
/*    * in identifying the data later on. */
/*    */
/*   void activate (const std::string& descriptor=""); */

/*   /** */
/*    * When the \p MeshData should be used, but was @e not activated */
/*    * prior to reading in a mesh, then the compatibility mode enables */
/*    * to still use this object as if the \p MeshData was active. */
/*    * The foreign node and element ids are simply assigned the */
/*    * indices used in \p libMesh.  Note that the compatibility mode */
/*    * should be used with caution, since the node and element */
/*    * indices in \p libMesh may be renumbered any time.  This */
/*    * \p MeshData always employs the current node and element ids, */
/*    * it does @e not create an image of ids when compatibility */
/*    * mode was activated. */
/*    */ 
/*   void enable_compatibility_mode (const std::string& descriptor=""); */

/*   /** */
/*    * Clears the data fields, but leaves the id maps */
/*    * untouched.  Useful for clearing data for a new */
/*    * data file.  Use \p slim() to delete the maps. */
/*    */
/*   void clear (); */

/*   /** */
/*    * Once the data is properly read from file, the id  */
/*    * maps can safely be cleared.  However, if this object */
/*    * should remain able to @e write nodal or element oriented  */
/*    * data to file, this method should better @e not be used. */
/*    * Use the appropriate \p bool to select the id map that */
/*    * should be cleared.  By default, both id maps are deleted. */
/*    */ 
/*   void slim (const bool node_id_map = true, */
/* 	     const bool elem_id_map = true); */

/*   /** */
/*    * Translates the @e nodal data contained in this object */
/*    * to \p data_values and \p data_names.  These two */
/*    * vectors are particularly suitable for use with */
/*    * the \p MeshBase::write method that takes nodal */
/*    * data.  E.g., the export method may be used for */
/*    * inspecting boundary conditions.  A reference */
/*    * to the mesh for which the data should be written */
/*    * has to be provided.  Note that this mesh @e has  */
/*    * to contain the nodes for which this \p MeshData  */
/*    * holds data.  I.e., \p out_mesh may only refer to  */
/*    * the \p MeshBase itself (that this \p MeshData belongs  */
/*    * to), or its \p BoundaryMesh, cf. \p Mesh.   */
/*    */ 
/*   void translate (const MeshBase& out_mesh, */
/* 		  std::vector<Number>& data_values, */
/* 		  std::vector<std::string>& data_names) const; */







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






/* protected: */


/*   //---------------------------------------------------------- */
/*   // read/write   Methods */
/*   /** 


/*  /** */
/*    * Read nodal/element oriented data in UNV format, */
/*    * either from an ASCII file or from a gzip'ed ASCII  */
/*    * file, using the C++ wrapper \p gzstream to \p zlib.h. */
/*    */ 
/*   //  void read_unv (const std::string& file_name); */
/*   void read_unv_elements (const std::string& file_name); */


/*   /** */
/*    * Actual implementation of reading nodal/element  */
/*    * oriented data in UNV format.  This has to be */
/*    * decoupled from \p read_unv() in order to allow */
/*    * reading both \p .unv and \p .unv.gz files. */
/*    */ 
/*   //  void read_unv_implementation (std::istream& in_file); */

/*   void read_unv_implementation_elements (std::istream& in_file); */






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
