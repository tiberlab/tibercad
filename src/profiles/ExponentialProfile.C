// $Id$

#include "tibercad/profiles/ExponentialProfile.h"
#include "tibercad/base/InitFailedException.h"

#include <limits>

using namespace std;


ExponentialProfile::ExponentialProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _max(0.0),
  _lambda(1),
  _type(onesided),
  _direction(1, 0, 0),
  _origin(0)
{
  _max = get_option("max", _max);
  get_option("origin", _origin);
  get_option("direction", _direction);
  _lambda = get_option("lambda", _lambda);
  double decay = 0.0;
  decay = get_option("decay_length", decay);
  if (decay > 0.0)
    _lambda = 1.0 / decay;

  _direction /= _direction.norm();

  string type = get_option("type", "onesided");
  if (type == "symmetric")
    _type = symmetric;
  else if (type == "continued")
    _type = continued;

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
ExponentialProfile::get_data(const Elem* , const Point& p) const
{
  double data = 0.0;

  // shift to origin
  Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  if ((_type == symmetric) || (xcoord >= 0))
    data = _max * exp(-fabs(xcoord) * _lambda);

  if ((_type == continued) && (xcoord < 0))
    data = _max;

  return data;
}
