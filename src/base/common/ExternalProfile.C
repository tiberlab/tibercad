// $Id$

#include "ExternalProfile.h"

#include "ExtProfile1D.h"
#include "UniformRandomAlloy.h"

using namespace std;


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
  else if (options.get_name() == "uniform_random")
    pr = new UniformRandomAlloy(options);

  return pr;
}


double
ExternalProfile::get_data(const Elem* elem) const
{
  double conc = 0.0;
  int n_nodes = elem->n_nodes();

  for (int i = 0; i < n_nodes; ++i)
    conc += get_data(elem, elem->point(i));


  return(conc / n_nodes);
}
