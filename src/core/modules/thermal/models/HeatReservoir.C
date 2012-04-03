 // $Id: HeatReservoir.C 2362 2011-02-21 01:02:31Z gromano $

#include "HeatReservoir.h"

#include "TiberModule.h"


void
HeatReservoir::do_init(void)
{
  get_parameter("temperature", _temperature);
}

void
HeatReservoir::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  set_coefficients(1.0, 0.0, _temperature);
}

