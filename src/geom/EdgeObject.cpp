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
 * \file EdgeObject.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/geom/EdgeObject.h"
#include "tibercad/io/Database.h"

#include <cassert>


EdgeObject*
EdgeObject::create(const ModelOptions& options)
{
  // In this case we cannot associate a default database
  Database db("", options.get_option("datafile", ""));
  db.set_section("");

  EdgeObject* eo = new EdgeObject(options);

  assert(eo != NULL);

  eo->set_database(db);

  return eo;
}



