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
 * \file EdgeObject.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_EDGEOBJECT_H
#define TC_EDGEOBJECT_H


#include "tibercad/physics/PhysicalObject.h"

class Material;


//! Description of an edge, i.e. (n-2)-D object.
class EdgeObject : public PhysicalObject
{

  public:

    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    ~EdgeObject(void) {};


    //! Create an edge object
    /*!
     * \param options options for this boundary
     */
    static EdgeObject* create(const ModelOptions& options) TBDLLOCAL;


  protected:

    //! Construct an edge object
    EdgeObject(const ModelOptions& options) : PhysicalObject(EDGE, options) {};


    // \copydoc PhysicalObject::do_init()
    //void do_init(void);


  private:

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------



#endif /* _EDGEOBJECT_H_ */
