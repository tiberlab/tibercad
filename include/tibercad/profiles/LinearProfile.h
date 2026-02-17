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
 * \file LinearProfile.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_LINEARPROFILE_H
#define TC_LINEARPROFILE_H


#include "tibercad/profiles/ExternalProfile.h"

#include "point.h"
#include "elem.h"

#include <vector>

//! Class to create a linear profile
class LinearProfile : public ExternalProfile
{
  public:

    //! Constructor
    LinearProfile(const ModelOptions& options);

    //! Destructor
    virtual ~LinearProfile(void);

    virtual double get_data(const libMesh::Elem* elem, const libMesh::Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

  private:

    double _min;
    double _max;

    double _distance;

    enum Type {
      onesided,
      symmetric,
      continued,
    };

    Type _type;

    libMesh::Point _direction;
    libMesh::Point _origin;

};


#endif //_LINEARPROFILE_H_
