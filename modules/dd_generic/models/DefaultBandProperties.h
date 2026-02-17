/*  
 * This file is part of the tiberCAD module dd_generic.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file DefaultBandProperties.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_DEFAULTBANDPROPERTIES_H
#define TC_DEFAULTBANDPROPERTIES_H

#include "CarrierProperties.h"


//! Base class for band parameter models
class DefaultBandProperties : public CarrierProperties
{

  public:

    //! Destructor
    virtual ~DefaultBandProperties(void) = default; 

    //! Creator method
    static DefaultBandProperties* create(const ModelOptions& options);


  protected:

    //! Constructor
    DefaultBandProperties(const ModelOptions& options);

  private:

};


#endif // TC_DEFAULTBANDPROPERTIES_H
