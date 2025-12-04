/*  
 * This file is part of the tiberCAD module wateringress.
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
 * \file WIBoundaryModel.h
 * \brief tiberCAD wateringress module header.
 *
 * \note This file is part of module wateringress.
 */


#ifndef _WIBOUNDARYMODEL_H_
#define _WIBOUNDARYMODEL_H_

#include "tibercad/physics/PhysicalModel.h"


class Elem;
class Point;


//! The base class for WI boundary conditions
class WIBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~WIBoundaryModel(void) {};

    //! Creator function
    static WIBoundaryModel* create(const MaterialBoundary* boundary,
        const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

    void get_coefficients(double& a, double& b, double& c);


  protected:

    //! Constructor
    WIBoundaryModel(const ModelOptions& options);

    void set_coefficients(double a, double b, double c);


  private:

    double _alpha;
    double _beta;
    double _gamma;

};



inline
WIBoundaryModel::WIBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _alpha(0),
  _beta(1),
  _gamma(0)
{
}



inline
void
WIBoundaryModel::get_coefficients(double& a, double& b, double& c)
{
  a = _alpha;
  b = _beta;
  c = _gamma;
}


inline
void
WIBoundaryModel::set_coefficients(double a, double b, double c)
{
  _alpha = a;
  _beta = b;
  _gamma = c;
}


#endif // _WIBOUNDARYMODEL_H_
