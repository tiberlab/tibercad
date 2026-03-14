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
 * \file ISE_Element_0D.h
 * \brief Private tiberCAD header.
 *
 * \internal
 */

#ifndef ISE_ELEMENT_0D_H_
#define ISE_ELEMENT_0D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Vertex.h"
#include <algorithm>


using namespace  std;


//! 0D Class Element.
/*!
  Point Definition.
*/
class TC_DLLOCAL ISE_Element_0D : public ISE_Element
{
 public:
  //! Constructor.
  /*!
    Assigns a Vertex pointer to the element.
  */
  ISE_Element_0D(ISE_Vertex* vertex_0);

  //! Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Element_0D();



 private:

  /*!
    Vertex Pointer.
  */
  ISE_Vertex* node;

};

inline
ISE_Element_0D::ISE_Element_0D(ISE_Vertex* vertex_0) :ISE_Element()
{
  node = vertex_0;
  element_nodes_id.resize(1);
  element_nodes_id[0] = node->get_node_id();

}



#endif /*ISE_ELEMENT_0D_H_*/
