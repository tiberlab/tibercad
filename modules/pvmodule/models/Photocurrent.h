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
 * \file Photocurrent.h
 * \brief tiberCAD pvmodule module header.
 *
 * \note This file is part of module pvmodule.
 */

#ifndef TC_PHOTOCURRENT_H
#define TC_PHOTOCURRENT_H

#include "tibercad/physics/PhysicalModel.h"

/*!
 * \brief Base class for photocurrent models
 *
 * The model has to implement a method that
 * returns the photocurrent density in function
 * of a given coordinate.
 */
class Photocurrent : public PhysicalModel
{

  public:

    virtual ~Photocurrent(void) {};

    double get_photocurrent(const libMesh::Elem* elem,
                            const libMesh::Point& p) const
    {
      return do_get_photocurrent(elem, p);
    };


  protected:

    //! Constructor
    Photocurrent(const ModelOptions& options)
      : PhysicalModel(options) {};


    virtual double do_get_photocurrent(const libMesh::Elem* elem,
                                       const libMesh::Point& p) const = 0;

  private:

};


#endif // TC_PHOTOCURRENT_H
