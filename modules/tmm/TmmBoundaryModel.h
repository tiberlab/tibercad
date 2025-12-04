/*  
 * This file is part of the tiberCAD module tmm.
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
 * \file TmmBoundaryModel.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */


#ifndef _TMMBOUNDARYMODEL_H_
#define _TMMBOUNDARYMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
#include "Tmm.h"


namespace libMesh
{
  class Elem;
  class Point;
}

//! The base class for thermal balance boundary conditions
class TmmBoundaryModel : public  PhysicalModel
{

public:

  //! Destructor
  ~TmmBoundaryModel(void) {};

  //! Creator function
  static TmmBoundaryModel* create(const MaterialBoundary* boundary,
    const ModelOptions& options);

  virtual void Calculate_M_Matrix(void) = 0;
  void set_elements(double, double, double, double);
  virtual double get_element(int);
  virtual double get_kr(void);
  virtual double get_steps(void);
  void set_dipole_elements(double, double);
  std::string read_type(void);








protected:

  //! Constructor
  TmmBoundaryModel(const ModelOptions& options);
  void write_type(std::string);





private:
  std::string typer;
  double _mmm00;
  double _mmm01;
  double _mmm10;
  double _mmm11;
  double __kr;
  double __steps;



};


inline
TmmBoundaryModel::TmmBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options)
{
}



#endif // _TMMBOUNDARYMODEL_H_
