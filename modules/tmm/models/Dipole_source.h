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
 * \file Dipole_source.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */

/*
 * incidentwave.h
 *
 *  Created on: 4 Oct 2021
 *      Author: pamiri
 */

#ifndef SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_
#define SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_

#include "TmmBoundaryModel.h"
#include "Tmm.h"

namespace libMesh
{
  class Elem;
}
class TC_DLLOCAL Dipole_source : public TmmBoundaryModel
{

  public:

    //! Constructor
    Dipole_source(const ModelOptions& options);

    //! Destructor
    ~Dipole_source(void) {};

    //! Creator function

    static Dipole_source* create(const ModelOptions& options);
    virtual void Calculate_M_Matrix(void);


  protected:

    virtual void do_init(void);




  private:
    double _kr;
    double _steps;


};



inline
Dipole_source::Dipole_source(const ModelOptions& options) :
  TmmBoundaryModel(options)
{
}



inline
Dipole_source*
Dipole_source::create(const ModelOptions& options)
{
  return new Dipole_source(options);
}

#endif /* SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_ */
