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
 * \file TightBindingModelInterface.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef _TIGHTBINDINGMODELINTERFACE_H_
#define _TIGHTBINDINGMODELINTERFACE_H_

#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/physics/PhysicalModel.h"

//!Class that contains objects necessary for Atomistic Tight Binding Calculations
class TightBindingModelInterface : public PhysicalModel
{
public:

  //!Constructor
  TightBindingModelInterface(const ModelOptions& options);

  //!Destructor
  ~TightBindingModelInterface(){};

  //!Create a new object
  //static TightBindingModelInterface* create(const std::string& name,  const ModelOptions& options = ModelOptions());


protected:


  virtual PhysicalModel* create_new (void) const = 0;

  virtual void do_init(void) = 0;

};


inline
TightBindingModelInterface::TightBindingModelInterface(const ModelOptions& options) :
  PhysicalModel(options)
{
}


//inline
//TightBindingModelInterface* TightBindingModelInterface::create(const std::string& name,  const ModelOptions& options)
//{
//  return dynamic_cast<TightBindingModelInterface*> ( PhysicalModel::create(name,options) );
//}


#endif
