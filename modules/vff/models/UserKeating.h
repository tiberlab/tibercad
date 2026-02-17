/*  
 * This file is part of the tiberCAD module vff.
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
 * \file UserKeating.h
 * \brief tiberCAD vff module header.
 *
 * \note This file is part of module vff.
 */

#ifndef TC_USERKEATING_H
#define TC_USERKEATING_H

#include "tibercad/base/tiber_dll.h"
#include "Keating.h"

//! User defined Keating model parameters
class TBDLLOCAL UserKeating : public Keating
{
public:

  //! Destructor
  ~UserKeating(void) {};

  //! Creator function
  static UserKeating* create(const ModelOptions& options);

  double get_alpha();

  //! Assign value to parameters
  void do_init(void);

protected:

private:

  UserKeating(const ModelOptions& options);

  void assign_alpha(void);

  void assign_beta(void);

  //! Parse parent material parameters. To be used if the models belongs to an alloy
  void assign_alpha_parents(void);

  void assign_beta_parents(void);


};

inline
UserKeating*
UserKeating::create(const ModelOptions& options)
{
  return new UserKeating(options);
}



#endif
