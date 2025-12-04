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
 * \file SingleBand.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */


#include "tibercad/physics/PhysicalModel.h"
#include "HamiltonianModel.h"
#include "tensor_value.h"
#include "vector_value.h"

class SingleBand : public HamiltonianModel 
{

  public:

    virtual ~SingleBand(void) {};

    static SingleBand* create(const ModelOptions& options);

  protected:

    SingleBand(const ModelOptions& options);

    void read_database(void);

    void do_init(void);

    void set_invmass_tensor(void);

};

inline
SingleBand*
SingleBand::create(const ModelOptions& options)
{
  return new SingleBand(options);
}


