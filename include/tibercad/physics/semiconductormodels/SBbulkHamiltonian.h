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
 * \file SBbulkHamiltonian.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _SBBULKHAMILTONIAN_H_
#define _SBBULKHAMILTONIAN_H_

#include "tibercad/physics/semiconductormodels/EFAbulkHamiltonian.h"
#include "tibercad/math/Tensor2.h"

#include <complex>
#include <vector>



//! A class that builds single band Hamiltonian
class SBbulkHamiltonian : public EFAbulkHamiltonian
{

 public:

  //! destructor
  ~SBbulkHamiltonian(void);

 

 

  virtual void calculate_Hamiltonian_k_par(void);

 

  virtual void calculate_Hamiltonian_gen(void); 


  //! Applies ONLY potential. Strain in applied only in derived classes
  virtual void apply_strain_and_potential(const Tensor2& strain_crystal, double el_potential);

   





 

 protected:

  //! default constructor
  SBbulkHamiltonian(const ModelOptions& options);

  //!band edge
  double edge;

  //!\f$ \left(\frac{1}{m}\right)_{ij} \f$ tensor
  Tensor2 imass;

  //!The only matrix element of the Hamiltonian
  MatrixElement single_band_ham;

 
  //! a pointer to a semiconductor that contains parameters
  //DDsemiconductor* semiconductor;


  virtual PhysicalModel* create_new(void) const = 0;

  virtual void do_init(void) = 0;

  virtual void do_print_info(void);

 private:

  

  
 
};


#endif
