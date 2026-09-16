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
 * \file ExternalProfile.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_EXTERNALPROFILE_H
#define TC_EXTERNALPROFILE_H

#include "tibercad/module/SolutionProvider.h"

#include "tibercad/module/TiberModelObject.h"
#include "tibercad/base/libMeshDefs.h"

class Elem;

/*!
 * \brief Base class for reading external data profiles
 */
class ExternalProfile : protected TiberModelObject
{
  public:

    //! Destructor
    virtual ~ExternalProfile(void);

    //! The creation method
    static ExternalProfile* create(const ModelOptions& options);

    //! Get the data at an element
    virtual double get_data(const Elem* elem) const;

    //! Get the data at a coordinate
    virtual double get_data(const Elem* elem, const Point& p) const;

    //! Get extremal values
    virtual std::pair<double, double> get_min_max(void) const;


  protected:

    //! Constructor
    ExternalProfile(const ModelOptions& options);


  private:

    //! Setup, called only if now derived class is instantiated
    void setup(void);

    //! The data provider
    SolutionProvider _data_source;

};



#endif //_EXTERNALPROFILE_H_
