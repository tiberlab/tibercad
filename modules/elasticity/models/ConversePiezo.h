/*  
 * This file is part of the tiberCAD module elasticity.
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
 * \file ConversePiezo.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef TC_CONVERSEPIEZO_H
#define TC_CONVERSEPIEZO_H

#include "BodyForceModel.h"

#include "tibercad/module/SimulationInterface.h"
#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"

class Elem;



//! The base class for Poisson boundary conditions
class TC_DLLOCAL ConversePiezo : public BodyForceModel
{

  public:
 
  //! Destructor
  ~ConversePiezo(void){};
  
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point) override;

  protected:

    //! Constructor
    ConversePiezo(const ModelOptions& options);

    //! Initialize
    virtual void do_init(void) override;


    /* This is not used here: */
     virtual void read_database(void) override;





  private:
  
  SimulationInterface* _simul;

  double _e33;
  double _e31;
  double _e15;
  
  ID ElFieldID;

  
};





#endif // TC_GRAYMODEL_H
