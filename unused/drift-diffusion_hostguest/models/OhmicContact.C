// $Id: OhmicContact.C 3542 2013-03-01 09:31:59Z maufder $

#include "OhmicContact.h"
#include "DriftDiffusionProperties.h"

#include "TiberModule.h"


OhmicContact::OhmicContact(const ModelOptions& options)
 : ElectricalContact(options)
{
  set_type(0, DIRICHLET);
  //set_type(0, NEUMANN);

  // as default we apply the potential
  set_type(1, DIRICHLET);
  set_type(2, DIRICHLET);

}



void
OhmicContact::do_compute(void)
{
  set_contact_fermilevel(get_bulk_dd_properties()->get_equilibrium_fermi_level());

  ElectricalContact::do_compute();
}
