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
 * \file EFAbulkModel.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _EFABULKMODEL_H_
#define _EFABULKMODEL_H_

#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/physics/semiconductormodels/EFAbulkHamiltonian.h"

//! A container class that contains a dynamical object of bulk Hamiltonian
class EFAbulkModel: public PhysicalModel
{

 public: 

  //!Destructor
  ~EFAbulkModel() = default;

  //!creates a new object
  static EFAbulkModel* create(const ModelOptions& options);

  EFAbulkHamiltonian* get_Hamiltonian_model(void) const;

 protected:

  //!Constructor
  EFAbulkModel(const ModelOptions& options);

  virtual PhysicalModel* create_new (void) const;


  virtual void do_init_alloy (const PhysicalModel *comp_A, const PhysicalModel *comp_B, double xa);

  virtual void do_init(void);
  
  virtual void do_print_info(void);

  virtual void prepare_submodels(void);



 private:

  //! a pointer to the bulk Hamiltonian
  EFAbulkHamiltonian*  _bulkHamiltonian = nullptr;


};


inline EFAbulkModel*  EFAbulkModel::create(const ModelOptions& options)
{
  return new EFAbulkModel(options);
}

inline  PhysicalModel* EFAbulkModel::create_new (void) const
{
  return new EFAbulkModel(get_options());
}

inline  EFAbulkHamiltonian* EFAbulkModel::get_Hamiltonian_model(void) const
{
  return  _bulkHamiltonian;
}

#endif
