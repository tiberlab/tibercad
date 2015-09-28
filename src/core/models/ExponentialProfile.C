// $Id$

#include "ExponentialProfile.h"
#include "InitFailedException.h"

#include <limits>

using namespace std;


ExponentialProfile::ExponentialProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _max(0.0),
  _lambda(1),
  _direction(1, 0, 0),
  _origin(0)
{
  _max = get_option("max", _max);
  get_option("origin", _origin);
  get_option("direction", _direction);
  _lambda = get_option("lambda", _lambda);
  _direction /= _direction.size();
}

ExponentialProfile::~ExponentialProfile(void)
{
}





pair<double, double>
ExponentialProfile::get_min_max(void) const
{
  return(make_pair(0.0, _max));
}




double
ExponentialProfile::get_data(const Elem* elem, const Point& p) const
{
  double data = _max;

  // shift to origin
  Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  if (xcoord > 0)
  {
    data = _max * exp(-xcoord * _lambda);
  }

  return data;
}
