/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file OhmicContact.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */


#ifndef TC_OHMICCONTACT_H
#define TC_OHMICCONTACT_H

#include "ElectricalContact.h"


class TBDLLOCAL OhmicContact : public ElectricalContact
{
  public:

    //! Create an ohmic contact
    static OhmicContact* create(const ModelOptions& options);

  protected:

    //! The constructor
    OhmicContact(const ModelOptions& options);


    //! Calculate all coefficients
    virtual void do_compute(void);


};


//
// inline
//

inline
OhmicContact*
OhmicContact::create(const ModelOptions& options)
{
  return new OhmicContact(options);
}






#endif // TC_OHMICCONTACT_H
