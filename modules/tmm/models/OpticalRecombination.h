/*  
 * This file is part of the tiberCAD module tmm.
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
 * \file OpticalRecombination.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */


#ifndef _OpticalRecombination_H_
#define _OpticalRecombination_H_

#include "Tmm.h"
#include "TmmDipoleSource.h"

namespace libMesh
{
  class Elem;
}

// Base class for InCoherence model
class TBDLEXPORT OpticalRecombination : public TmmDipoleSource
{

  public:

  virtual ~OpticalRecombination(void) {};
  static OpticalRecombination* create(const ModelOptions& options);
  virtual void calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda);

  
protected:
  
    OpticalRecombination(const ModelOptions& options);
	
	virtual void do_init(void){};





  private:

   double _recombination_rate;
   
       //! The generation model
    std::vector<SimulationInterface*> _recombination_model;

    //! The solution ID of the generation models variable
    std::vector<ID> _recombination_id;

};

inline
OpticalRecombination::OpticalRecombination(const ModelOptions& options) :
  TmmDipoleSource(options),
  _recombination_rate(0.0)
{
}

inline
OpticalRecombination*
OpticalRecombination::create(const ModelOptions& options)
{
  return new OpticalRecombination(options);
}



#endif // _POLARIZATIONMODEL_H_
