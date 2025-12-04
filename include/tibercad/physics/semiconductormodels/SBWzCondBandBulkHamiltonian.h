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
 * \file SBWzCondBandBulkHamiltonian.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _SBWzCondBandBulkHamiltonian_h_
#define _SBWzCondBandBulkHamiltonian_h_


#include "tibercad/physics/semiconductormodels/SBbulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/WzSemiconductor.h"
#include "tibercad/physics/semiconductormodels/SBCondBandBulkHamiltonian.h"

//! A class to calculate single band Hamiltonian of wurtzite material
class SBWzCondBandBulkHamiltonian : public SBCondBandBulkHamiltonian
{

 public:

  //!Destructor 
  ~SBWzCondBandBulkHamiltonian(){};


  
  virtual void apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential);

  
  static SBWzCondBandBulkHamiltonian* create(const ModelOptions& options);
  

 private:

 

  const WzSemiconductor::WzDDparameters*  wz_par; 

 protected:

  //!Constructor
  SBWzCondBandBulkHamiltonian(const ModelOptions& options)
    : SBCondBandBulkHamiltonian(options) {};

  virtual PhysicalModel* create_new(void) const ;

 
  virtual void calculate_for_init(void);


};

inline PhysicalModel* SBWzCondBandBulkHamiltonian::create_new() const
{
  return new SBWzCondBandBulkHamiltonian(get_options());
}

inline SBWzCondBandBulkHamiltonian* SBWzCondBandBulkHamiltonian::create(const ModelOptions& options)
{
  return new SBWzCondBandBulkHamiltonian(options);
}


#endif
