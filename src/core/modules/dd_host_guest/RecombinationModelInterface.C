#include "RecombinationModelInterface.h"

void
RecombinationModelInterface::set_coupling(const std::string& coupling)
{
  if (coupling != "")
    _coupling = coupling;
  else
    _coupling = get_option("coupling", "eh");

}
