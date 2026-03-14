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
 * \file FourierModel.h
 * \brief tiberCAD boltzmann module header.
 *
 * \note This file is part of module boltzmann.
 */


#ifndef TC_FOURIERMODEL_H
#define TC_FOURIERMODEL_H

#include "HeatTransportModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TC_DLLOCAL FourierModel : public HeatTransportModel
{

  public:
 
    //! Creator function
    static FourierModel* create(const ModelOptions& options);

   //! Destructor
    ~FourierModel(void) {};

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
  FourierModel(const ModelOptions& options);
  
};


inline
PhysicalModel*
FourierModel::create_new(void) const
{
  return new  FourierModel(get_options());
}

inline
FourierModel*
FourierModel::create(const ModelOptions& options)
{
  return new  FourierModel(options);
}
#endif // TC_GRAYMODEL_H
