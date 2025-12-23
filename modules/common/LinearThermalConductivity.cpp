/*  
 * This file is part of the tiberCAD module common.
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
 * \file LinearThermalConductivity.C
 * \brief tiberCAD common module implementation.
 *
 * \note This file is part of module common.
 */



#include "tibercad/physics/Material.h"
#include "LinearThermalConductivity.h"
#include "tibercad/base/libMeshDefs.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/base/ModelOptions.h"

#include "tibercad/module/TiberModule.h"



using namespace std;


LinearThermalConductivity::LinearThermalConductivity(const ModelOptions& options):ThermalConductivityModel(options)
{
   kx0 = 0.0;
   kz0 = 0.0;
   mx  = 0.0;
   mz  = 0.0;
   z0  = 0.0;
}



void
LinearThermalConductivity::do_init(void)
{

  const ModelOptions& options = get_options();


   kx0 = 0.0;
   kz0 = 0.0;
   mx  = 0.0;
   mz  = 0.0;
   z0  = 0.0;

 
   get_parameter("kx0",kx0);
   get_parameter("mx",mx);
   get_parameter("kz0",kz0);
   get_parameter("mz",mz);
   get_parameter("z0",z0);



  
}



void 
LinearThermalConductivity::calculate(const libMesh::Elem* , const libMesh::Point& point, double )
{
   
  double x = point(0);
  double y = point(1);
  double z = point(2);

  double kx = kx0  + mx * (z-z0);
  double kz = kz0  + mz * (z-z0);


  this->set_thermal_conductivity(libMesh::RealGradient(kx, kx, kz));
  //_kappa(0,0) = kx;
  //_kappa(1,1) = kx;
  //_kappa(2,2) = kz;


}
