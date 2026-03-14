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
 * \file ISE_Face.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef ISE_FACE_H_
#define ISE_FACE_H_

#include "ISE_Vertex.h"
#include "ISE_Edge.h"

#include <iostream>
#include <cstdlib>


using  namespace  std ;

//! Face Class
/*! Face Class. Contains edges and orientations.
 */
class TC_DLLOCAL ISE_Face
{
 public:

  //! Constructor.
  /*!
    Assigns edge pointers and orientations vector.
  */
  ISE_Face(vector <ISE_Edge*> fcs_edg, vector <bool> neg_edg) ;
	
 


  //! Gets edge belonging to a face.
  /*!
   *  Returns Edge pointer.
   *  e=0 -> first face edge; e=1 -> second face edge; ...
   *  An overflow control is included ('e' value limited).
   */

  ISE_Edge* get_edge (unsigned int e);

 
  /*!
   * Return the 'e' edge state.
   *  (true == negative edge; false == positive edge).
   */
  bool get_edge_state (unsigned int e);


  /*!
   * Returns face edges number.
   */
  unsigned int get_face_edges_size();

  //! Virtual Destructor.
  /*!
    Dummy. 
  */	
  virtual ~ISE_Face();
	
 private:

  /*!
    Face Edges.
  */
  vector <ISE_Edge*> face_edges;

  /*! 
    Orientation vector.
  */
  vector <bool> negative_face_edges;
};


inline ISE_Face::ISE_Face(vector <ISE_Edge*> fcs_edg, vector <bool> neg_edg)
{
	
  for (unsigned int e=0; e < fcs_edg.size(); e++)
  {
    face_edges.push_back(fcs_edg[e]);
    negative_face_edges.push_back(neg_edg[e]);
  }

}



inline ISE_Edge* ISE_Face::get_edge (unsigned int e)
{

  if (e < face_edges.size())
  {
    return face_edges[e];
  }
  else
  {
    cerr << "The 'e' value for the 'ISE_Edge* ISE_Face::get_edge(unsigned int e)' method is incorrect" << endl;
    std::exit(1);
  }
}


inline bool ISE_Face::get_edge_state (unsigned int e)
{

  if (e < negative_face_edges.size())
  {
    return negative_face_edges[e];
  }
  else
  {
    cerr << "The 'e' value for the 'bool ISE_Face::get_edge_state(unsigned int e)' method is incorrect" << endl;
    exit(1);
  }
}



inline unsigned int ISE_Face::get_face_edges_size()
{
  return face_edges.size();
}

#endif /*ISE_FACE_H_*/
