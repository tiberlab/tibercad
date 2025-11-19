 // $Id: HeatReservoir.C 2362 2011-02-21 01:02:31Z gromano $

#include "Diffusive.h"
#include "tibercad/base/SimulationOptions.h"

#include "tibercad/module/TiberModule.h"


void
Diffusive::do_init(void)
{
  //get_parameter("deltaT", _deltaT);

  get_parameter("p", _p);
 //
  //get_parameter("temperature", _temperature);

}

void
Diffusive::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{

  //set_deltaT(_deltaT);

  //set_periodicity(_periodicity);


  set_coefficients(1.0, 0.0, _p);

}

