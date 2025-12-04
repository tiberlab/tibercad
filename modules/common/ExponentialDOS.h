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
 * \file ExponentialDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef _EXPONENTIALDOS_H_
#define _EXPONENTIALDOS_H_


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT ExponentialDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~ExponentialDOS(void) {};


    //! Creator function
    static ExponentialDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    ExponentialDOS(const ModelOptions& options);

    virtual void do_init(void);


    //! Get occupied states and derivative
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double E, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;



  private:

    //! The tail parameter
    double _alpha;

    //! Calculates the value under the integral
    double _get_value(double e, double E, double kT) const;


};

//
// inline methods
//

inline
ExponentialDOS*
ExponentialDOS::create(const ModelOptions& options)
{
  return new ExponentialDOS(options);
}


#endif // _EXPONENTIALDOS_H_
