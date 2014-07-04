// $Id$

#include "ExternalProfile.h"

#include "ExtProfile1D.h"


ExternalProfile::ExternalProfile(const ModelOptions& options) :
  TiberModelObject(options)
{
}


ExternalProfile::~ExternalProfile(void)
{
}


ExternalProfile*
ExternalProfile::create(const ModelOptions& options)
{
  ExternalProfile* pr = nullptr;

  if (options.get_name() == "1d_profile")
    pr = new ExtProfile1D(options);

  return pr;
}
