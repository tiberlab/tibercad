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
 * \file ExtProfile1D.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _EXTPROFILE1D_H_
#define _EXTPROFILE1D_H_


#include "tibercad/profiles/ExternalProfile.h"

#include "point.h"
#include "elem.h"

#include <vector>

//! Class to read an external line profile
class ExtProfile1D : public ExternalProfile
{
  public:

    //! Constructor
    ExtProfile1D(const ModelOptions& options);

    //! Destructor
    virtual ~ExtProfile1D(void);

    virtual double get_data(const Elem* elem, const Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

  private:

    void _read_source(void);

    std::vector<double> _x_coord;
    std::vector<double> _values;
    std::vector<unsigned int> _addressing;

    double _min;
    double _max;

    Point _direction;
    Point _origin;
    double _scale;
    double _data_scale;

};


#endif //_EXTPROFILE1D_H_
