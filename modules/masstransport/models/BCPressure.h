/*  
 * This file is part of the tiberCAD module masstransport.
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
 * \brief tiberCAD masstransport module header.
 *
 * \note This file is part of module masstransport.
 */


#ifndef TC_BCPRESSURE_H
#define TC_BCPRESSURE_H

#include "MTBoundaryModel.h"
#include <string>
using namespace std;




//! The base class for Poisson boundary conditions
class BCPressure : public MTBoundaryModel
{

  public:

    //! Destructor
    ~BCPressure(void) {};

    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point);


  protected:

    //! Constructor
    BCPressure(const ModelOptions& options);

    //! Initialize
    virtual void do_init(void);


  private:

    //! The relative pressure
    double _relative_pressure = 0; 

    //! The molecule to be considered
    string _mass_tran_spec = "H2O"; // this can be defined in the input file

};



inline
BCPressure::BCPressure(const ModelOptions& options) :
  MTBoundaryModel(options)
{
}



#endif // TC_BCPRESSURE_H

