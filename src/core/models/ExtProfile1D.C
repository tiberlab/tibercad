// $Id$

#include "ExtProfile1D.h"
#include "InitFailedException.h"

#include <limits>
#include <fstream>

using namespace std;


ExtProfile1D::ExtProfile1D(const ModelOptions& options) :
  ExternalProfile(options),
  _min(numeric_limits<double>::max()),
  _max(numeric_limits<double>::min()),
  _direction(1, 0, 0),
  _origin(0),
  _scale(1.0),
  _data_scale(1.0)
{
  get_option("origin", _origin);
  get_option("direction", _direction);
  _direction /= _direction.norm();
  _scale = get_option("length_scaling", _scale);
  _data_scale = get_option("data_scaling", _data_scale);

  _max = get_option("max_value", _max);
  _min = get_option("min_value", _min);

  _read_source();

}

ExtProfile1D::~ExtProfile1D(void)
{
}


void
ExtProfile1D::_read_source(void)
{

  // a line buffer
  const size_t buf_len = 256;
  char buf[buf_len];

  string file = get_option("source", "");
  ifstream is(file);
  if (is.fail() || !is.good())
    throw InitFailedException("Cannot read 1D profile "
        "from \"" + file + "\"");


  size_t i = 0;
  while (is.good())
  {
    if (i == _x_coord.size())
    {
      size_t n_new = _x_coord.size() + 100;
      _x_coord.reserve(n_new);
      _values.reserve(n_new);
    }

    is.getline(buf, buf_len);
    if ((buf[0] != '#') && (buf[0] != '%') && (buf[0] != '/'))
    {
      istringstream in(buf);

      double l, s;
      if (in >> l)
      {
        if (in >> s)
        {
          _x_coord.push_back(_scale * l);
          double data = _data_scale * s;
          _values.push_back(data);

          if (data > _max)
            _max = data;
          if (data < _min)
            _min = data;

          i++;
        }
      }
    }
  }
  is.close();

  _x_coord.resize(_x_coord.size());
  _values.resize(_values.size());

  const int N = 50;
  double dx = (_x_coord.back() - _x_coord.front()) / N;
  _addressing.resize(N);

  unsigned int ctr = 0;
  for (unsigned int i = 0; i < N; ++i)
  {
    double x = _min + i * dx;

    while ((_x_coord[ctr] < x) && (ctr < _x_coord.size()))
      ++ctr;

    _addressing[i] = ctr;
  }
}


pair<double, double>
ExtProfile1D::get_min_max(void) const
{
  return(make_pair(_min, _max));
}




double
ExtProfile1D::get_data(const Elem* , const Point& p) const
{
  double data;

  // shift to origin
  Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  // calculate index for adressing array
  double dx = (_x_coord.back() - _x_coord.front()) / _addressing.size();
  int index = max(0, static_cast<int>(floor((xcoord - _x_coord.front()) / dx)));
  index = min(index, static_cast<int>(_addressing.size() - 1));

  unsigned int ctr = _addressing[index];

  while ((ctr < _x_coord.size()) && (_x_coord[ctr] < xcoord))
      ++ctr;

  if (ctr == 0)
    data = _values[0];
  else if (ctr == _x_coord.size())
    data = _values.back();
  else
  {
    double x1 = _x_coord[ctr - 1];
    double x2 = _x_coord[ctr];

    double frac = (xcoord - x1) / (x2 - x1);
    data = (_values[ctr] - _values[ctr - 1]) * frac + _values[ctr - 1];
  }

  return data;
}
