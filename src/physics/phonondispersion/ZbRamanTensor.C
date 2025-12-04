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
 * \file ZbRamanTensor.C
 * \brief Internal tiberCAD code.
 *
 * \internal
 */

#include "ZbRamanTensor.h"
#include "getpot.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "RotatedCrystal.h"  
#include "tibercad/module/SimulationEnvironment.h" 
//--------------------------------------------------------//

ZbRamanTensor::ZbRamanTensor(const ModelOptions& options)
 : RamanTensor(options)
{
}


void  ZbRamanTensor::read_database(void)
{
 
  const Material* mat = get_material();

  GetPot data((mat->get_database()).get_data_file());

  _raman_d = data("raman_d", 1.0);

 
}


void  ZbRamanTensor::do_init(void)
{

   const ModelOptions& options = get_options();

   _raman_d =  options.get_option("raman_d",_raman_d);

   _raman_tensor.resize(3);

   for (unsigned int n = 0;n<3;n++)
         {_raman_tensor[n]= 0 ;}   
        
       _raman_tensor[0](3,2) = _raman_d;
       _raman_tensor[1](3,1) = _raman_d;
       _raman_tensor[2](2,1) = _raman_d;

       Material* mat = get_material();

       const RotatedCrystal&   cr = mat->get_rotated_crystal();
     
       rotate_to_calculation_system(cr.RotMatrix);
 


}
