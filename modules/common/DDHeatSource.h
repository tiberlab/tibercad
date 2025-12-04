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
 * \file DDHeatSource.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef _DDHEATSOURCE_H_
#define _DDHEATSOURCE_H_

#include "tibercad/physics/misc/HeatSourceModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/base/tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL DDHeatSource : public HeatSourceModel
{

  public:
 
     //! Destructor
  ~DDHeatSource(void) {};
  
  //! Creator function
  static DDHeatSource* create(const ModelOptions& options);

  virtual void calculate(const Elem* elem, const Point& point);
    

  protected:

    //! Initialize
    virtual void do_init(void);

    /* In some cases it might be useful to reimplement this: */
    // virtual void do_init_interface(const PhysicalModel* comp_A,
    //         const PhysicalModel* comp_B);


    /* This is not used here: */
    // virtual void read_database(void);


    /* We do not use this here: */
    // virtual void read_interface_database(void);



  private:
  
    double _heat_source;

  enum heat_variables
    {
      EJOULE = 0,
      HJOULE,
      RECHEAT,
      EPELTH,
      HPELTH,
      WNX,
      WPX
    };
  
  //!Heat source variables for drift diffusion
  std::set<ID> ID_set;
  
  //!Variable map
  std::map<ID,ID> var_map;

  //!Pointer to drift diffusion simulation
  SimulationInterface* _simul;

  //! Constructor
    DDHeatSource(const ModelOptions& options);
  
};




inline
DDHeatSource*
DDHeatSource::create(const ModelOptions& options)
{
  return new  DDHeatSource(options);
}




#endif // _GRAYMODEL_H_
