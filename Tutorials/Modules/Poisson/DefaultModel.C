// $Id$

#include "DefaultModel.h"

TIBER_MODULE(DefaultModel, default)

using namespace std;

void
DefaultModel::do_init(void)
{
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

  set_charge_density(0);
}
