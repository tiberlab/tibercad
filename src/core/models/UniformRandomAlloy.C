// $Id$

#include "UniformRandomAlloy.h"
#include "InitFailedException.h"
#include "Messages.h"

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

