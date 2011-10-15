// $Id: FluxBoundary.C 2362 2011-02-21 01:02:31Z gromano $

#include "FluxBoundary.h"

//TIBER_MODULE(FluxBoundary,boltzmann_bnd, flux)

void
FluxBoundary::do_init(void)
{
  get_parameter("heat_flux", _heat_flux);
}

void
FluxBoundary::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(0, 1, _heat_flux);
}

