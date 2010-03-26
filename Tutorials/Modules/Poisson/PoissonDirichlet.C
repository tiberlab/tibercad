// $Id$

#include "PoissonDirichlet.h"



TIBER_MODULE(PoissonDirichlet, dirichlet)


void
PoissonDirichlet::do_init(void)
{
  get_parameter("potential", _potential);
}


void
PoissonDirichlet::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(1, 0, _potential);
}
