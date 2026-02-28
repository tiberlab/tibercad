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
 * \file TmmDipoleSource.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */


#ifndef _TmmDipoleSource_H_
#define _TmmDipoleSource_H_

#include "tibercad/physics/PhysicalModel.h"
#include "elem.h"
// #include "Tmm.h"

using namespace std;

//! This is the base class for the TMM bulk physical model
class TmmDipoleSource : public PhysicalModel
{

public:


  //! Destructor
  virtual ~TmmDipoleSource(void) {};

  const double& get_emission_power(void) const;


  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda) {};






protected:

  //! Constructor
  TmmDipoleSource(const ModelOptions& options);


  void set_emission_power(const double& emission_power);


private:


  double _emission_power;



};




#endif // TC_TMMBULKMODEL_H
