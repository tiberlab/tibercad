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
 * \file Pyropolarization.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_PYROPOLARIZATION_H
#define TC_PYROPOLARIZATION_H

#include "tibercad/physics/misc/PolarizationModel.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/tiber_dll.h"

// Base class for charge density models
class  TBDLLOCAL Pyropolarization: public PolarizationModel
{

  public:
  
   virtual ~Pyropolarization(void) {};
   
   static Pyropolarization* create(const ModelOptions& options);
  
  protected:

   Pyropolarization(const ModelOptions& options);

   virtual void do_init(void);

   virtual void read_database(void);

   virtual void do_calculate(const libMesh::Elem*, const libMesh::Point&) {};

  private:

    double _Pz;

    //! Initialize P from given Pz
    void _initP(void);

};


inline
Pyropolarization::Pyropolarization(const ModelOptions& options) :
  PolarizationModel(options),
  _Pz(0)
{
}


inline
Pyropolarization*
Pyropolarization::create(const ModelOptions& options)
{
  return new Pyropolarization(options);
}


#endif // TC_PIEZOPOLARIZATION_H
