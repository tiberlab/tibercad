/*  
 * This file is part of the tiberCAD module pvmodule.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file ElementaryCell.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */


#ifndef TC_ELEMENTARYCELL_H
#define TC_ELEMENTARYCELL_H

#include "tibercad/physics/PhysicalModel.h"

#include <ostream>


// Base class for elementary cell models
class ElementaryCell : public PhysicalModel
{

  public:

    virtual ~ElementaryCell(void) {};

    //! Write the netlist
    /*!
     * This method writes the elementary cell netlist into the
     * circuit file. Note that \c next_free must be a valid unused
     * node id at return.
     *
     * \param top_node the index of the top layer node
     * \param bottom_node the index of the bottom layer node
     * \param next_free the next free circuit node, will be updated during
     *        the call
     * \param area the cell area in cm^2
     * \param elem the mesh element pointer
     * \param p the coordinates of the elementary cell
     * \param os the stream to write to
     */
    void write_netlist(unsigned int top_node, unsigned int bottom_node,
                       unsigned int& next_free,
                       double area,
                       const libMesh::Elem* elem,
                       const libMesh::Point& p,
                       std::ostream& os) const
    { do_write_netlist(top_node, bottom_node, next_free, area, elem, p, os); };


  protected:

    ElementaryCell(const ModelOptions& options);


    //! Write the netlist
    /*!
     * This method writes the elementary cell netlist into the
     * circuit file.
     *
     * \param top_node the index of the top layer node
     * \param bottom_node the index of the bottom layer node
     * \param next_free the next free circuit node, will be updated during
     *        the call
     * \param area the cell area in cm^2
     * \param elem the mesh element pointer
     * \param p the coordinates of the elementary cell
     * \param os the stream to write to
     */

    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  std::ostream& os) const = 0;

  private:


};


inline
ElementaryCell::ElementaryCell(const ModelOptions& options) :
  PhysicalModel(options)
{
}



#endif // TC_ELEMENTARYCELL_H
