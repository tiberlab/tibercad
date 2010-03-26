// $Id$

#include "EdgeObject.h"
#include "Database.h"

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



