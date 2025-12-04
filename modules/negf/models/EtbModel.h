/*  
 * This file is part of the tiberCAD module negf.
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
 * \file EtbModel.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */


#include "HamiltonianModel.h"

class EtbModel : public HamiltonianModel 
{

  public:

    virtual ~EtbModel(void){};

    static EtbModel* create(const ModelOptions& options);

  protected:

    EtbModel(const ModelOptions& options);

    void do_init(void);


};

inline
EtbModel*
EtbModel::create(const ModelOptions& options)
{
  return new EtbModel(options);
}


