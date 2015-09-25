// $Id$

#include "LinearProfile.h"
#include "InitFailedException.h"

#include <limits>

using namespace std;


LinearProfile::LinearProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _min(0.0),
  _max(0.0),
  _direction(1, 0, 0),
  _origin(0)
{
  _min = get_option("min", _min);
  _max = get_option("max", _max);
  get_option("origin", _origin);
  get_option("direction", _direction);
  _distance = get_option("distance", _distance);
  _direction /= _direction.size();
}

LinearProfile::~LinearProfile(void)
{
}





pair<double, double>
LinearProfile::get_min_max(void) const
{
  return(make_pair(_min, _max));
}




double
LinearProfile::get_data(const Elem* elem, const Point& p) const
{
  double data = _max;

  // shift to origin
  Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  if (xcoord > 0)
  {
    if (xcoord < _distance)
      data = _max - xcoord / _distance * (_max - _min);
    else
      data = _min;
  }

  return data;
}
