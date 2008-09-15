// $Id$

#include "DSSCContact.h"


DSSCContact*
DSSCContact::create(const std::string& name,
    const ModelOptions& options)
{
  DSSCContact* ct = NULL;

  ct = new DSSCContact();

  if (ct != NULL)
    ct->set_options(options);

  if (name == "ohmic")
    ct->is_cathode() = true;

  return ct;
}



void
DSSCContact::do_init(void)
{
  std::string s(get_options().get_option("voltage", ""));
  set_potential(check_and_register(s, 0.0));
}



void
DSSCContact::calculate_current(void)
{
}
