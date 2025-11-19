// $Id$

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

