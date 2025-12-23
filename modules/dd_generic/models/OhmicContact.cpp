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
 * \file OhmicContact.C
 * \brief tiberCAD dd_generic module implementation.
 *
 * \note This file is part of module dd_generic.
 */


#include "OhmicContact.h"
#include "DriftDiffusionProperties.h"

#include "tibercad/module/TiberModule.h"


OhmicContact::OhmicContact(const ModelOptions& options)
 : ElectricalContact(options)
{
  

}

void
OhmicContact::do_init(void)
{
  for (unsigned int i = 0; i <= n_known_carriers(); i++)
    set_type(i, DIRICHLET);

  ElectricalContact::do_init();
}

void
OhmicContact::do_compute(void)
{
  set_contact_fermilevel(get_bulk_dd_properties()->get_equilibrium_fermi_level());

  ElectricalContact::do_compute();
}
