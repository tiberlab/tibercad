// $Id$

#include "PoissonDirichlet.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_dirichlet.so
TIBER_MODULE(PoissonDirichlet, contact, dirichlet)


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
