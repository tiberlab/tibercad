// $Id: Clamp.C 2124 2010-10-22 14:00:17Z gromano $

#include "BuiltInStrain.h"

#include "TiberModule.h"




void
BuiltInStrain::do_init(void)
{
  
  RealTensor strain(0);
  get_parameter("strain", strain);

  
  set_strain_source(strain);

}


