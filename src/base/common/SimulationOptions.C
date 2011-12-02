// $Id$


#include "SimulationOptions.h"
#include "ModelOptions.h"
#include "Variable.h"





double
SimulationOptions::temperature = 300.0;

double&
SimulationOptions::temp = temperature;

double&
SimulationOptions::T = temperature;

bool
SimulationOptions::incomplete_ionization = false;

int
SimulationOptions::_verbose = 1;



int
SimulationOptions::verbose(void)
{
  return _verbose;
}


void
SimulationOptions::initialize(const ModelOptions& opts)
{
  std::string temp = opts.get_option("temperature", "300.0");
  Variable::check_and_register(temp, temperature);

  incomplete_ionization = opts.get_option("incomplete_ionization", true);

  _verbose = opts.get_option("verbose", 2);
}

