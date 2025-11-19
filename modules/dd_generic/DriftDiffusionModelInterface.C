// $Id: DriftDiffusionModelInterface.C 2882 2011-07-18 10:09:45Z maufder $

#include "DriftDiffusionModelInterface.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/physics/Material.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/base/RuntimeException.h"

const double
DriftDiffusionModelInterface::T0 = Constants::k_B * 300.0;


DriftDiffusionProperties&
DriftDiffusionModelInterface::get_bulk_driftdiffusionproperties(void) const
{
  if (get_material() == NULL)
    throw RuntimeException("Drift-Diffusion model \'" + get_name()
      + "\' has no bulk material assigned but needs it.");

  return dynamic_cast<DriftDiffusionProperties&>(*
      get_material()->get_model(get_simulator_id()));

}

DriftDiffusionProperties&
DriftDiffusionModelInterface::get_driftdiffusionproperties(void) const
{
  if (get_owner() == NULL)
    throw RuntimeException("Drift-Diffusion model \'" + get_name()
      + "\' has no bulk material assigned but needs it.");

  return dynamic_cast<DriftDiffusionProperties&>(*
      get_owner()->get_model(get_simulator_id()));

}
