// $Id: PoissonNeumann.C 2029 2010-07-21 13:01:58Z maufder $

#include "PoissonNeumann.h"

#include "TiberModule.h"



void
PoissonNeumann::do_init(void)
{
  // taken in C/m^2
  get_parameter("displacement", _field);
}


void
PoissonNeumann::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(0, 1, _field);
}
