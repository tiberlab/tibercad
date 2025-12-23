/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file CompositeProfile.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/profiles/CompositeProfile.h"
#include "tibercad/base/InitFailedException.h"

#include <limits>

using namespace std;


CompositeProfile::CompositeProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _peak(1.0),
  _offset(0.0)
{
  _peak = get_option("peak_value", _peak);
  _offset = get_option("offset", _offset);

  ModelOptions& opts = get_options();
  auto it = opts.submodels_begin("profile");
  for ( ; it != opts.submodels_end("profile"); ++it)
  {
    if (!it->second.find_option("max"))
      it->second.set_option("max", 1.0);
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

  return(_peak * data + _offset);
}
