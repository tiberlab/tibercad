/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file SBuserHamiltonian.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _SBUSERHAMILTONIAN_H_
#define _SBUSERHAMILTONIAN_H_
//! A class that builds single band Hamiltonian

#include <complex>
#include <vector>
#include "tibercad/physics/semiconductormodels/SBbulkHamiltonian.h"
#include "tibercad/physics/PhysicalModel.h"


class SBuserHamiltonian : public SBbulkHamiltonian
{

 public:

  //! destructor
  ~SBuserHamiltonian(void);


  static  SBuserHamiltonian* create(const ModelOptions& options);

  

 protected:

  //! default constructor
  SBuserHamiltonian(const ModelOptions& options);



  virtual PhysicalModel* create_new(void) const;

 

  virtual void do_init(void);

 
  virtual void do_init_alloy (const PhysicalModel *comp_A, const PhysicalModel *comp_B, double xa) {};

  


 private:

 

  
 
};

inline PhysicalModel* SBuserHamiltonian::create_new() const
{
  return new SBuserHamiltonian(get_options());
}





inline SBuserHamiltonian* SBuserHamiltonian::create(const ModelOptions& options)
{
  return new SBuserHamiltonian(options);
}
#endif
