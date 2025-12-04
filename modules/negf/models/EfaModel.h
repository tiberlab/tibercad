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
 * \file EfaModel.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */


#include "HamiltonianModel.h"

class EfaModel : public HamiltonianModel 
{

  public:

    virtual ~EfaModel(void){};

    static EfaModel* create(const ModelOptions& options);

  protected:

    EfaModel(const ModelOptions& options);

    void do_init(void);

};

inline
EfaModel*
EfaModel::create(const ModelOptions& options)
{
  return new EfaModel(options);
}


