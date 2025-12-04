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
 * \file GrayModel.h
 * \brief tiberCAD boltzmann module header.
 *
 * \note This file is part of module boltzmann.
 */


#ifndef _GRAYMODEL_H_
#define _GRAYMODEL_H_

#include "HeatTransportModel.h"

#include "point.h"
#include "tensor_value.h"
#include "vector_value.h"
#include "tibercad/base/tiber_dll.h"

class Elem;


//! The base class for Poisson boundary conditions
class TBDLLOCAL GrayModel : public HeatTransportModel
{

  public:

 //! Creator function
  static GrayModel* create(const ModelOptions& options);

  //! Destructor
  ~GrayModel(void) {};


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
    GrayModel(const ModelOptions& options);

 //  struct options
//   {
    
//     double max_iter;
//     double max_error; //!< Max tollerance for self-consistent loop
//     vector<ID> spec;
//     std::vector<double> d_omega;
//     std::vector<Point> directions;
//     std::vector<Point> dir;
//     std::vector<double> theta_vec;
//     std::vector<double> phi_vec;
//     ID theta_slices;
//     ID phi_slices;
   
//   };

//  options myopts;


};



inline
PhysicalModel*
GrayModel::create_new(void) const
{
  return new  GrayModel(get_options());
}

inline
GrayModel*
GrayModel::create(const ModelOptions& options)
{
  return new  GrayModel(options);
}




#endif // _GRAYMODEL_H_
