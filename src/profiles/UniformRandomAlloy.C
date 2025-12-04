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
 * \file UniformRandomAlloy.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/profiles/UniformRandomAlloy.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/io/Messages.h"

#include "elem.h"

using namespace std;

UniformRandomAlloy::UniformRandomAlloy(const ModelOptions& options) :
  ExternalProfile(options),
  _min(0.0),
  _max(1.0),
  _rnd_generator(0)
{

  _min = get_option("min_value", _min);
  _max = get_option("max_value", _max);

  if (has_option("site_density"))
    _site_density = get_option("site_density", _site_density);
  else
    throw InitFailedException("Uniform random alloy profile needs the density "
        "of available sites (in mesh units) as option.");

  if (has_option("mean_composition"))
    _mean_composition = get_option("mean_composition", _mean_composition);
  else
    throw InitFailedException("Uniform random alloy profile needs the "
        "mean composition as option.");

  int seed = get_option("random_generator_seed",
      static_cast<int>(time(NULL) * random_device()()));

  std::ostringstream os;
  os << "Initializing  MT19937 random generator with seed " << seed;
  Messages::info(os.str());

  _rnd_generator.seed(seed);
}

UniformRandomAlloy::~UniformRandomAlloy(void)
{
}


pair<double, double>
UniformRandomAlloy::get_min_max(void) const
{
  return(make_pair(_min, _max));
}


double
UniformRandomAlloy::get_data(const Elem* elem) const
{
  double volume = elem->volume();

  int n_sites = _site_density * volume;
  binomial_distribution<int> binomial(n_sites, _mean_composition);
  double data = binomial(_rnd_generator);
  data /= n_sites;

  if (data > _max) data = _max;
  else if (data < _min) data = _min;

  return(data);
}


double
UniformRandomAlloy::get_data(const Elem* elem, const Point& p) const
{
  return(get_data(elem));
}

