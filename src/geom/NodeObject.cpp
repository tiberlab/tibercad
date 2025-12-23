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
 * \file NodeObject.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/geom/NodeObject.h"
#include "tibercad/io/Database.h"

#include <cassert>


NodeObject*
NodeObject::create(const ModelOptions& options)
{
  // In this case we cannot associate a default database
  Database db("", options.get_option("datafile", ""));
  db.set_section("");

  NodeObject* no = new NodeObject(options);

  assert(no != NULL);

  no->set_database(db);

  return no;
}

