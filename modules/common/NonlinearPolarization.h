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
 * \file NonlinearPolarization.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef _PIEZOPOLARIZATION_H_
#define _PIEZOPOLARIZATION_H_

#include "tibercad/physics/misc/PolarizationModel.h"
#include "tibercad/physics/StrainInterface.h"
#include "tibercad/base/tiber_dll.h"

class Elem;
class Point;

/*!
 * \brief Nonlinear polarization model including spontaneous polarization
 *
 * based on Phys. Rev. B 88, 121304(R) (2013)
 */
class TBDLLOCAL NonlinearPolarization: public PolarizationModel
{

  public:

    virtual ~NonlinearPolarization(void) {};

    static NonlinearPolarization* create(const ModelOptions& options);


  protected:

    NonlinearPolarization(const ModelOptions& options);

    virtual void do_init(void);

    virtual void read_database(void);

    virtual void do_calculate(const Elem* elem, const Point& point);

  private:

    //! The strain simulation
    StrainInterface _strain;

    //! The spontaneous polarization (along [0 0 0 1])
    double _Psp;

    //! Piezoelectric modulus \f$e_{33}\f$ (wurtzite)
    double _e33;

    //! Piezoelectric modulus \f$e_{31}\f$ (wurtzite)
    double _e31;

    //! Piezoelectric modulus \f$e_{15}\f$ (wurtzite) or \f$e_{14}\f$ (zincblende)
    union
    {
        double _e15;
        double _e14;
    };

    /*!
     * \brief 2nd order coefficients
     *
     * The coefficients in \c _coeff are, according to PRB 88:
     * 2a 2b c 2d 2e f g h
     */
    std::vector<double> _2nd_order_coeff;
};


inline
NonlinearPolarization::NonlinearPolarization(const ModelOptions& options) :
  PolarizationModel(options),
  _Psp(0),
  _e33(0),
  _e31(0),
  _e15(0)
{
  _2nd_order_coeff = {0, 0, 0, 0, 0, 0, 0, 0};
}


inline
NonlinearPolarization*
NonlinearPolarization::create(const ModelOptions& options)
{
  return new NonlinearPolarization(options);
}


#endif // _PIEZOPOLARIZATION_H_
