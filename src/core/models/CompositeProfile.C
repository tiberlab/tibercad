// $Id$

#include "CompositeProfile.h"
#include "InitFailedException.h"

#include <limits>

using namespace std;


CompositeProfile::CompositeProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _peak(1.0)
{
  _peak = get_option("peak_value", _peak);

  ModelOptions& opts = get_options();
  auto it = opts.submodels_begin("profile");
  for ( ; it != opts.submodels_end("profile"); ++it)
  {
    if (!it->second.find_option("max"))
      it->second.set_option("max", 1.0);
    it->second.print_all();
    _profiles.push_back(ExternalProfile::create(it->second));
  }
  
}

CompositeProfile::~CompositeProfile(void)
{
  for (auto&& p : _profiles)
    delete p;
}





pair<double, double>
CompositeProfile::get_min_max(void) const
{
  double min = 1.0, max = 1.0;
  for (auto&& p : _profiles)
  {
    auto minmax = p->get_min_max();
    min *= minmax.first;
    max *= minmax.second;
  }

  return(make_pair(min, max));
}




double
CompositeProfile::get_data(const Elem* elem, const Point& p) const
{
  double data = 1.0;

  for (auto&& it : _profiles)
    data *= it->get_data(elem, p);

  return(_peak * data);
}
