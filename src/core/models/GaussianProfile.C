// $Id$

#include "GaussianProfile.h"
#include "InitFailedException.h"

#include <limits>

using namespace std;


GaussianProfile::GaussianProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _max(0.0),
  _sigma(1.0),
  _direction(1, 0, 0),
  _origin(0)
{
  _max = get_option("max", _max);
  get_option("origin", _origin);
  get_option("direction", _direction);
  _sigma = get_option("sigma", _sigma);
  _direction /= _direction.size();
}

GaussianProfile::~GaussianProfile(void)
{
}





pair<double, double>
GaussianProfile::get_min_max(void) const
{
  return(make_pair(0.0, _max));
}




double
GaussianProfile::get_data(const Elem* elem, const Point& p) const
{
  double data = _max;

  // shift to origin
  Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  if (xcoord > 0)
  {
    data = _max * exp(-xcoord * xcoord / (_sigma * _sigma));
  }

  return data;
}
