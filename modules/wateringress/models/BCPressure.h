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
 * \file BCPressure.h
 * \brief tiberCAD wateringress module header.
 *
 * \note This file is part of module wateringress.
 */


#ifndef _BCPRESSURE_H_
#define _BCPRESSURE_H_

#include "WIBoundaryModel.h"




//! The base class for Poisson boundary conditions
class BCPressure : public WIBoundaryModel
{

  public:

    //! Destructor
    ~BCPressure(void) {};

    //! Creator function
    static BCPressure* create(const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point);


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

    //! Constructor
    BCPressure(const ModelOptions& options);

    //! The relative humidity
    double _relative_humidity = 0;

};



inline
BCPressure::BCPressure(const ModelOptions& options) :
  WIBoundaryModel(options)
{
}



inline
BCPressure*
BCPressure::create(const ModelOptions& options)
{
  return new BCPressure(options);
}



#endif // _BCPRESSURE_H_
