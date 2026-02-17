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
 * \file IncoherentModel.h
 * \brief tiberCAD tmm module header.
 *
 * \note This file is part of module tmm.
 */


#ifndef TC_INCOHERENCE_H
#define TC_INCOHERENCE_H

#include "Tmm.h"
#include "TmmBulkModel.h"

namespace libMesh
{
  class Elem;
}

// Base class for InCoherence model
class TBDLEXPORT IncoherentModel : public TmmBulkModel
{

  public:

  virtual ~IncoherentModel(void) {};
  
  const double& get_Incoherent_Index(void) const;
  static IncoherentModel* create(const ModelOptions& options);

  
protected:
  
    IncoherentModel(const ModelOptions& options);

    void set_Incoherent_Index(const double& Incoheret_Index);




  private:

   double _Incoheret_Index;

};




#endif // TC_POLARIZATIONMODEL_H
