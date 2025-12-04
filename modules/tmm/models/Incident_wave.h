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
 * \file Incident_wave.h
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
class TBDLLOCAL Incident_wave : public TmmBoundaryModel
{

  public:

    //! Destructor
    ~Incident_wave(void) {};

    //! Creator function


    static Incident_wave* create(const ModelOptions& options);
    virtual void Calculate_M_Matrix(void);


  protected:

    virtual void do_init(void);




  private:

    //! Constructor
    Incident_wave(const ModelOptions& options);


};



inline
Incident_wave::Incident_wave(const ModelOptions& options) :
  TmmBoundaryModel(options)
{
}



inline
Incident_wave*
Incident_wave::create(const ModelOptions& options)
{
  return new Incident_wave(options);
}

#endif /* SRC_CORE_MODULES_TMM_MODELS_INCIDENTWAVE_H_ */
