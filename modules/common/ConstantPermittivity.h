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
 * \file ConstantPermittivity.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef _CONSTANTPERMITTIVITY_H_
#define _CONSTANTPERMITTIVITY_H_

#include "tibercad/physics/misc/PermittivityModel.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/io/Database.h"


// Base class for charge density models
class  TBDLLOCAL ConstantPermittivity: public PermittivityModel
{

  public:

    virtual ~ConstantPermittivity(void) {};

    static ConstantPermittivity* create(const ModelOptions& options);


  protected:

    ConstantPermittivity(const ModelOptions& options);

    virtual void do_init(void);

    virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point){};

    virtual void read_database(void);

  private:

    libMesh::RealVectorValue _permittivity_diag;

};


inline
ConstantPermittivity::ConstantPermittivity(const ModelOptions& options) :
  PermittivityModel(options),
  _permittivity_diag(0)
{
}


inline
ConstantPermittivity*
ConstantPermittivity::create(const ModelOptions& options)
{
  return new ConstantPermittivity(options);
}


#endif // _PIEZOPOLARIZATION_H_
