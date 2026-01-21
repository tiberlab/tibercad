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
 * \file LinearProfile.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/profiles/LinearProfile.h"
#include "tibercad/base/InitFailedException.h"

#include <limits>

using namespace std;


LinearProfile::LinearProfile(const ModelOptions& options) :
  ExternalProfile(options),
  _min(0.0),
  _max(0.0),
  _distance(1.0),
  _type(onesided),
  _direction(1, 0, 0),
  _origin(0)
{
  _min = get_option("min", _min);
  _max = get_option("max", _max);
  get_option("origin", _origin);
  get_option("direction", _direction);
  _distance = get_option("distance", _distance);
  _direction /= _direction.norm();

  string type = get_option("type", "onesided");
  if (type == "symmetric")
    _type = symmetric;
  else if (type == "continued")
    _type = continued;

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
LinearProfile::get_data(const libMesh::Elem*, const libMesh::Point& p) const
{
  double data = 0.0;

  // shift to origin
  libMesh::Point point(p - _origin);

  // project onto _direction
  double xcoord = point * _direction;

  if (std::abs(xcoord) < _distance)
    data = _max - xcoord / _distance * (_max - _min);
  else
    data = _min;

  if (xcoord < 0)
    if (_type == continued)
      data = _max;
    else if (_type == onesided)
      data = 0.0;

  return data;
}
