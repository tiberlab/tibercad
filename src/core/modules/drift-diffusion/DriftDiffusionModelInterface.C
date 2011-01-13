// $Id$

#include "DriftDiffusionModelInterface.h"
#include "DriftDiffusionProperties.h"
#include "PhysicalObject.h"
#include "Constants.h"

const double
DriftDiffusionModelInterface::T0 = Constants::k_B * 300.0;


DriftDiffusionProperties&
DriftDiffusionModelInterface::get_driftdiffusionproperties(void)
{
  return dynamic_cast<DriftDiffusionProperties&>(*
      get_owner()->get_model(get_simulator_id()));
}
