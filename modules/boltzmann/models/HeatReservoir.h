/*  
 * This file is part of the tiberCAD module boltzmann.
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
 * \file HeatReservoir.h
 * \brief tiberCAD boltzmann module header.
 *
 * \note This file is part of module boltzmann.
 */


#ifndef TC_HEATRESERVOIR_H
#define TC_HEATRESERVOIR_H

#include "BoltzmannBoundaryModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"



//! The base class for Poisson boundary conditions
class TC_DLLOCAL HeatReservoir : public BoltzmannBoundaryModel
{

  public:

    //! Destructor
    ~HeatReservoir(void) {};

    //! Creator function
    static HeatReservoir* create(const ModelOptions& options);

   //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point);

    //! Calculate for a point on the given side
    //virtual void get_periodicity(const Point& point);

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


    //! Create a new object of the same type
    virtual PhysicalModel* create_new(void) const;


  private:

    //! Constructor
    HeatReservoir(const ModelOptions& options);

  double _temperature;
};



inline
HeatReservoir::HeatReservoir(const ModelOptions& options) :
  BoltzmannBoundaryModel(options),
  _temperature(0)
{
}



inline
HeatReservoir*
HeatReservoir::create(const ModelOptions& options)
{
  return new HeatReservoir(options);
}



inline
PhysicalModel*
HeatReservoir::create_new(void) const
{
  return new HeatReservoir(get_options());
}

#endif // TC_POISSONDIRICHLET_H
