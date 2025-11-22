#include "ZbFreeDynamicalMatrix.h"
#include "getpot.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "RotatedCrystal.h"  
//--------------------------------------------------------//
void  ZbFreeDynamicalMatrix::read_database(void)
{
 
  const Material* mat = get_material();

  GetPot data((mat->get_database()).get_data_file());

  w0 = data("w0", 0.0);
 
}

//---------------------------------------------------------//



void  ZbFreeDynamicalMatrix::do_init(void)
{

   const ModelOptions& options = get_options();

   w0 = options.get_option("w0",w0);
  
   double value = w0 * w0 /(8065.6 * 8065.6);
   _dynamical_matrix(1,1) = value;                                 
   _dynamical_matrix(2,2) = value;
   _dynamical_matrix(3,3) = value;
  
}

void  ZbFreeDynamicalMatrix::re_init(void)
{
}
