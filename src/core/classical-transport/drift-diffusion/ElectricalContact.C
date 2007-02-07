// $Id$

#include "OhmicContact.h"
#include "SchottkyContact.h"
#include "FermiLevelPinning.h"

#include "ElectricalContact.h"


ElectricalContact*
ElectricalContact::create(const std::string& name,
    const ModelOptions& options)
{
  ElectricalContact* ct = NULL;

  if (name == "ohmic")
    ct = OhmicContact::create();
  else if (name == "schottky")
    ct = SchottkyContact::create();
  else if (name == "pinning")
    ct = FermiLevelPinning::create();
  
  if (ct != NULL)
    ct->set_options(options);

  return ct;
}


void
ElectricalContact::do_init(void)
{
  if (get_options().get_option("zero_field", false))
    set_zero_derivative_bc(DriftDiffusionDefs::POTENTIAL);
  if (get_options().get_option("zero_grad_fermi_e", false))
    set_zero_derivative_bc(DriftDiffusionDefs::FERMIE);
  if (get_options().get_option("zero_grad_fermi_h", false))
    set_zero_derivative_bc(DriftDiffusionDefs::FERMIH);

  if (get_options().find_option("value"))
    set_simulation_voltage(get_options().get_option("value", 0.0));
}


