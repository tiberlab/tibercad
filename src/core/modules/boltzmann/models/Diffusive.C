 // $Id: HeatReservoir.C 2362 2011-02-21 01:02:31Z gromano $

#include "Diffusive.h"
#include "SimulationOptions.h"

TIBER_MODULE(Diffusive,boltzmann_bnd,diffusive)

void
Diffusive::do_init(void)
{
  //get_parameter("deltaT", _deltaT);

  //get_parameter("periodicity", _periodicity);

  //get_parameter("temperature", _temperature);

}

void
Diffusive::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{

  //set_deltaT(_deltaT);

  //set_periodicity(_periodicity);

  //set_coefficients(1.0, 0.0, SimulationOptions::temperature);


}

