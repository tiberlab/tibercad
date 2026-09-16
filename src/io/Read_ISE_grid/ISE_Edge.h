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
 * \file ISE_Edge.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef ISE_EDGE_H_
#define ISE_EDGE_H_
#include "ISE_Vertex.h"

using  namespace  std ;

//! Edge Class. 
/*!
  Contains two Vertex pointers.
*/
class TC_DLLOCAL ISE_Edge
{
 public:
  //!  Constructor 
  /*! 
    Assigns two Vertex Pointers.
  */
  ISE_Edge(ISE_Vertex*  id_1,  ISE_Vertex*   id_2) ;


  //!  Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Edge();
	
  //! Get node belonging to an  edge
  /*! Returns  a  pointer to  \c ISE_Vertex  belonging to an Edge. 
   *  v=1 -> first edge vertex; v=2 -> second edge vertex;
   */ 
  ISE_Vertex*  get_vertex(unsigned int  v); 
	
 



 private:

  /*!
    First node pointer.
  */
  ISE_Vertex*  node_1;

  /*!
    Second node pointer. 
  */
  ISE_Vertex*  node_2;

};




inline
ISE_Edge::ISE_Edge( ISE_Vertex*  id_1, ISE_Vertex*    id_2)
{
  node_1 = id_1;
  node_2 = id_2;
}

inline ISE_Vertex* 
ISE_Edge::get_vertex(unsigned int  v) 
{
  if  (v==1)
  {return node_1;}
  else
  {return node_2;}
}



#endif /*ISE_EDGE_H_*/
