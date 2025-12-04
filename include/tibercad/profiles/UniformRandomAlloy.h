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
 * \file UniformRandomAlloy.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _UNIFORMRANDOMALLOY_H_
#define _UNIFORMRANDOMALLOY_H_


#include "tibercad/profiles/ExternalProfile.h"

#include <random>


//! Class to create a uniform random distribution
class UniformRandomAlloy : public ExternalProfile
{
  public:

    //! Constructor
    UniformRandomAlloy(const ModelOptions& options);

    //! Destructor
    virtual ~UniformRandomAlloy(void);

    virtual double get_data(const Elem* elem) const override;

    virtual double get_data(const Elem* elem, const Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

  private:


    double _min;
    double _max;

    double _mean_composition;

    //! The available site density in mesh units
    double _site_density;

    //! Random number generator engine
    mutable std::mt19937 _rnd_generator;

};


#endif // _UNIFORMRANDOMALLOY_H_
