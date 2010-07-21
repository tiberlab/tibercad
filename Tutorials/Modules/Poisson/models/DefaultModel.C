// $Id$

#include "DefaultModel.h"
#include "MyChargeDensityModel.h"

// The first string is the class name, the second one
// is the type of the model (here it is a bulk model),
// the third one is the specific model implementation.
// The library name will then be bulk_default.so
TIBER_MODULE(DefaultModel, bulk, default)

using namespace std;

void
DefaultModel::do_init(void)
{
  SubmodelIterator it = submodels_begin("charge_density");
  if (it != submodels_end("charge_density"))
    _charge_density = dynamic_cast<ChargeDensityModel*>(it->second);
}


void
DefaultModel::read_database(void)
{
}


void
DefaultModel::do_calculate(void)
{
  RealTensor& eps = get_permittivity();
  eps(0,0) = eps(1,1) = eps(2,2) = 1.0;

  get_polarization().zero();

  double charge = 0.0;
  if (_charge_density != NULL)
    charge = _charge_density->get_charge_density(get_element(), get_point());

  set_charge_density(charge);
}
