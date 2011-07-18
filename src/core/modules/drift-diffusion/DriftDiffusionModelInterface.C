// $Id$

#include "DriftDiffusionModelInterface.h"
#include "DriftDiffusionProperties.h"
#include "Material.h"
#include "Constants.h"
#include "RuntimeException.h"

const double
DriftDiffusionModelInterface::T0 = Constants::k_B * 300.0;


DriftDiffusionProperties&
DriftDiffusionModelInterface::get_driftdiffusionproperties(void) const
{
  if (get_material() == NULL)
    throw RuntimeException("Drift-Diffusion model \'" + get_name()
      + "\' has no bulk material assigned but needs it.");

  return dynamic_cast<DriftDiffusionProperties&>(*
      get_material()->get_model(get_simulator_id()));
}
