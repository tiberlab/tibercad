// $Id: OhmicContact.C 3542 2013-03-01 09:31:59Z maufder $

#include "OhmicContact.h"
#include "DriftDiffusionProperties.h"

#include "TiberModule.h"


OhmicContact::OhmicContact(const ModelOptions& options)
 : ElectricalContact(options)
{
  

}

void
OhmicContact::do_init(void)
{
  ElectricalContact::do_init();

  for (unsigned int i = 0; i < n_known_carriers(); i++)
    set_type(i, DIRICHLET);
}

void
OhmicContact::do_compute(void)
{
  //set_contact_fermilevel(get_bulk_dd_properties()->get_equilibrium_fermi_level());

  ElectricalContact::do_compute();
}
