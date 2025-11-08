// $Id: PoissonDirichlet.C 3864 2014-06-13 10:24:26Z maufder $

#include "PoissonDirichlet.h"

#include "TiberModule.h"



void
PoissonDirichlet::do_init(void)
{
  get_parameter("potential", _potential);
}


void
PoissonDirichlet::calculate(const Elem* elem, unsigned int side,const Point& point)
{
  set_coefficients(1, 0, _potential);
}
