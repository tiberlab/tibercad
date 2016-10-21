// $Id$

#include "ConstantBodyForce.h"
#include "Material.h"

#include "TiberModule.h"


using namespace std;


ConstantBodyForce::ConstantBodyForce(const ModelOptions& options):BodyForceModel(options)
{
}

void
ConstantBodyForce::do_init(void)
{
  libMesh::RealGradient force_source(0);
  get_parameter("F", force_source);
  set_force_source(force_source);

  libMesh::RealTensor dummy_tens(0);
  set_strain_source(dummy_tens);
  set_stress_source(dummy_tens);
}




 
