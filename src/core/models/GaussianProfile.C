// $Id$

#include "GaussianProfile.h"
#include "InitFailedException.h"

#include <limits>

using namespace std;


GaussianProfile::GaussianProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _max(0.0),
  _sigma(1.0),
  _type(onesided),
  _direction(1, 0, 0),
  _origin(0)
{
  _max = get_option("max", _max);
  get_option("origin", _origin);
  get_option("direction", _direction);
  _sigma = get_option("sigma", _sigma);
  _direction /= _direction.size();
  
  string type = get_option("type", "onesided");
  if (type == "symmetric")
    _type = symmetric;
  else if (type == "continued")
    _type = continued;
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
  double data = 0.0;

  // shift to origin
  Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  if ((_type == symmetric) || (xcoord >= 0))
    data = _max * exp(-xcoord * xcoord / (_sigma * _sigma));

  if ((_type == continued) && (xcoord < 0))
    data = _max;

  return data;
}
