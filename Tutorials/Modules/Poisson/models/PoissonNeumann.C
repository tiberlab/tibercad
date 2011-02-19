// $Id: PoissonNeumann.C 2029 2010-07-21 13:01:58Z maufder $

#include "PoissonNeumann.h"



// The first string is the class name, the second one
// is the type of the model (here it is a contact model),
// the third one is the specific model implementation.
// The library name will then be contact_neumann.so
TIBER_MODULE(PoissonNeumann, contact, neumann)


void
PoissonNeumann::do_init(void)
{
  get_parameter("normal_field", _field);
}


void
PoissonNeumann::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(0, 1, _field);
}
