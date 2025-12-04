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
 * \file SBZbCondBandBulkHamiltonian.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _SBZbCondBandBulkHamiltonian_h_
#define _SBZbCondBandBulkHamiltonian_h_


#include "tibercad/physics/semiconductormodels/ZbSemiconductor.h"
#include "tibercad/physics/semiconductormodels/SBCondBandBulkHamiltonian.h"
#include "tibercad/physics/PhysicalModel.h"

//! A class to calculate single band Hamiltonian of zinc-blende material
class SBZbCondBandBulkHamiltonian : public SBCondBandBulkHamiltonian
{

 public:

  //!Destructor
  ~SBZbCondBandBulkHamiltonian(){};
  
  //! Gamma or X or L
  std::string min_name;

  //! minima number 
  /*!
    O for Gamma; 
    0,1,2 for X;
    0,1,2,3 for L.
  */ 
  unsigned int min_number;

  
  virtual void apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential);

  static SBZbCondBandBulkHamiltonian* create(const ModelOptions& options);

 private:

  const ZbSemiconductor::ZbDDparameters* zb_par;

 protected:

  //!Constructor
  SBZbCondBandBulkHamiltonian(const ModelOptions& options)
    : SBCondBandBulkHamiltonian(options) {};

  virtual PhysicalModel* create_new(void) const;


  virtual void calculate_for_init(void);

  virtual void do_print_info(void);


};

inline PhysicalModel* SBZbCondBandBulkHamiltonian::create_new() const
{
  return new SBZbCondBandBulkHamiltonian(get_options());
}

inline SBZbCondBandBulkHamiltonian* SBZbCondBandBulkHamiltonian::create(const ModelOptions& options)
{
  return new SBZbCondBandBulkHamiltonian(options);
}


#endif
