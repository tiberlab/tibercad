// $Id: LinearThermalConductivity.C 2457 2011-03-06 23:52:12Z gromano $


#include "Material.h"
#include "LinearThermalConductivity.h"
#include "libMeshDefs.h"
#include "PhysicalModel.h"
#include "ModelOptions.h"

#include "TiberModule.h"



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
