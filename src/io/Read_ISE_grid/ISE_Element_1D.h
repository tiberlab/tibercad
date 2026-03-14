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
 * \file ISE_Element_1D.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef ISE_ELEMENT_1D_H_
#define ISE_ELEMENT_1D_H_

#include "ISE_Element.h"
#include "ISE_Vertex.h"


//! 1D Element Class.
/*!
  Segment Definition.
*/
class TC_DLLOCAL ISE_Element_1D: public  ISE_Element
{
 public:

  //!  Constructor
  /*!
    Assigns two vertex pointers to the element.
  */
  ISE_Element_1D(ISE_Vertex*  vertex_0,  ISE_Vertex*  vertex_1);

  //!Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Element_1D();



 private:

  /*!
    Vertex Pointer.
  */
  ISE_Vertex*  node_1;

  /*!
    Vertex Pointer.
  */
  ISE_Vertex*  node_2;

};

inline
ISE_Element_1D::ISE_Element_1D(ISE_Vertex*  vertex_0, ISE_Vertex*  vertex_1 ):ISE_Element()
{
  node_1 = vertex_0;
  node_2 = vertex_1;

  element_nodes_id.resize(2);
  element_nodes_id[0] = node_1->get_node_id();
  element_nodes_id[1] = node_2->get_node_id();

}



#endif /*ISE_ELEMENT_1D_H_*/
