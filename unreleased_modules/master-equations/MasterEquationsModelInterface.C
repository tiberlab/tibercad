

#include "MasterEquationsModelInterface.h"
#include "MasterEquationsProperties.h"
#include "Material.h"
#include "Constants.h"
#include "RuntimeException.h"

const double
MasterEquationsModelInterface::T0 = Constants::k_B * 300.0;


MasterEquationsProperties&
MasterEquationsModelInterface::get_masterequationsproperties(void) const
{
  if (get_material() == NULL)
    throw RuntimeException("Master-Equations model \'" + get_name()
      + "\' has no bulk material assigned but needs it.");

  return dynamic_cast<MasterEquationsProperties&>(*
      get_material()->get_model(get_simulator_id()));
}
