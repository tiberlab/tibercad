 // $Id: HeatReservoir.C 2362 2011-02-21 01:02:31Z gromano $

#include "Periodic.h"
#include "SimulationOptions.h"

#include "TiberModule.h"


void
Periodic::do_init(void)
{
  get_parameter("deltaT", _deltaT);

  get_parameter("periodicity", _periodicity);

  //get_parameter("temperature", _temperature);

}

void
Periodic::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{

  set_deltaT(_deltaT);

  set_periodicity(_periodicity);

  set_coefficients(1.0, 0.0, SimulationOptions::temperature);


}

