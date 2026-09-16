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
 * \file SBCondBandBulkHamiltonian.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */



#ifndef TC_SBCONDBULKHAMILTONIAN_H
#define TC_SBCONDBULKHAMILTONIAN_H

#include "tibercad/physics/semiconductormodels/SBbulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/Semiconductor.h"


class Semiconductor;

//! A clas that builds a single band Hamiltonian for a conduction band of a crystal
class SBCondBandBulkHamiltonian: public SBbulkHamiltonian
{
 public:


 ~SBCondBandBulkHamiltonian(void);

  
 virtual void apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential) = 0;
  
  

 virtual void set_temperature(double Temperature);


 protected:

  SBCondBandBulkHamiltonian(const ModelOptions& options);

 virtual PhysicalModel* create_new(void) const = 0 ;


 virtual void do_init(void);

 virtual void prepare_submodels(void);

  

 //! a pointer to a semiconductor that contains parameters
 Semiconductor* semiconductor;

 //!calculate everything we need from semiconductor model
 virtual void calculate_for_init(void) = 0;

 private:


 

};
//----------------------------------------------------------------------------------------//

inline SBCondBandBulkHamiltonian::SBCondBandBulkHamiltonian(const ModelOptions& options)
  : SBbulkHamiltonian(options),
    semiconductor(NULL)
{
}



inline 
void SBCondBandBulkHamiltonian::set_temperature(double Temperature)
{
  semiconductor->set_temperature(Temperature);
  semiconductor->apply_temperature();
  calculate_for_init();
}


#endif
