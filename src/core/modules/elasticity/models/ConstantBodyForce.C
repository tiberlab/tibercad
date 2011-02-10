// $Id$

#include "ConstantBodyForce.h"
#include "Material.h"


// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so

TIBER_MODULE(ConstantBodyForce, BodyForce, constant)

using namespace std;


ConstantBodyForce::ConstantBodyForce(const ModelOptions& options):BodyForceModel(options)
{
}

void
ConstantBodyForce::do_init(void)
{
  RealGradient force_source(0);
  get_parameter("F", force_source);
  set_force_source(force_source);

  RealTensor dummy_tens(0);
  set_strain_source(dummy_tens);
  set_stress_source(dummy_tens);
}




 
