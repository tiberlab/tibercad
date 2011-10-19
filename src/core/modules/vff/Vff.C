#include "Vff.h"
#include "Messages.h"

TIBER_MODULE(Vff, MODULE_NAME)


Vff*
Vff::_this = NULL;


Vff::Options::Options(void)
: boundary_conditions("all_around"),
  substrate_plane("z")
{
}

Vff::Vff(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
}


Vff::~Vff(void)
{
  // there's nothing to be done
}


Vff*
Vff::create(const ModelOptions& options)
{
  return new Vff(options);
}
