/*  
 * This file is part of the tiberCAD module common.
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
 * \file Piezopolarization.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_PIEZOPOLARIZATION_H
#define TC_PIEZOPOLARIZATION_H

#include "tibercad/physics/misc/PolarizationModel.h"
#include "tibercad/module/SolutionProvider.h"
#include "tibercad/base/tiber_dll.h"


// Basic Piezopolarization model
class TC_DLLOCAL Piezopolarization: public PolarizationModel
{

  public:

    virtual ~Piezopolarization(void) {};

    static Piezopolarization* create(const ModelOptions& options);


  protected:

    Piezopolarization(const ModelOptions& options);

    virtual void do_init(void);

    virtual void read_database(void);

    virtual void do_calculate(const libMesh::Elem* elem, const libMesh::Point& point);

  private:

    //! The strain simulation
    SolutionProvider _strain;

    //! Piezoelectric modulus \f$e_{33}\f$ (wurtzite)
    double _e33;

    //! Piezoelectric modulus \f$e_{31}\f$ (wurtzite)
    double _e31;

    //! Piezoelectric modulus \f$e_{15}\f$ (wurtzite) or \f$e_{14}\f$ (zincblende)
    union
    {
        double _e15;
        double _e14;
    };
};


inline
Piezopolarization::Piezopolarization(const ModelOptions& options) :
  PolarizationModel(options),
  _e33(0),
  _e31(0),
  _e15(0)
{
}


inline
Piezopolarization*
Piezopolarization::create(const ModelOptions& options)
{
  return new Piezopolarization(options);
}


#endif // TC_PIEZOPOLARIZATION_H
