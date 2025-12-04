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
 * \file HeatSourceModel.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _HEATSOURCEMODEL_H_
#define _HEATSOURCEMODEL_H_

#include "tibercad/physics/PhysicalModel.h"

#include "tensor_value.h"

//#undef  TIBER_MODULE_PREFIX
//#define TIBER_MODULE_PREFIX thermal_conductivity

//typedef double Real;
#include "tibercad/base/libMeshDefs.h"
//class Elem;
//class Point;

using namespace std;

//! The base class for Poisson boundary conditions
class HeatSourceModel : public PhysicalModel
{

  public:

    //! Destructor
    ~HeatSourceModel(void) {};

     //! Creator function
   static HeatSourceModel* create(const ModelOptions& options);

  Real get_heat_source(void) const;

  virtual void calculate(const Elem* elem, const Point& point){};

  protected:

    //! Constructor
  HeatSourceModel(const ModelOptions& options);

  void set_heat_source(Real heat_source);

  private:

  Real _heat_source;

};

inline
Real
HeatSourceModel::get_heat_source(void) const
{

  return _heat_source;
}


inline 
void 
HeatSourceModel::set_heat_source(Real heat_source)
{
  _heat_source = heat_source;
}



inline
HeatSourceModel::HeatSourceModel(const ModelOptions& options) :
PhysicalModel(options)
{
}




#endif // _THERMALCONDUCTIVITYMODEL_H_
