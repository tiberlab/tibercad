// $Id$

#include "tibercad/physics/PhysicalObject.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/io/Database.h"

#include <cassert>


PhysicalObject::PhysicalObject(ObjectType type, const ModelOptions& options)
  : TiberModelObject(options),
    _type(type),
    _is_initialized(false)
{
  _database = new Database();
}


PhysicalObject::~PhysicalObject(void)
{
  ModelMap::iterator it(_models.begin());
  const ModelMap::iterator end(_models.end());
  for ( ; it != end; ++it)
    destroy(it->second);

  _models.clear();

  delete _database;
}



void
PhysicalObject::add_model(PhysicalModel* model, ID simulator_id)
{
  assert(simulator_id != 0);

  if (model != NULL)
  {
    ModelMap::iterator it = _models.find(simulator_id);
    if (it != _models.end())
    {
      destroy(it->second);
      it->second = model;
    }
    else
      _models[simulator_id] = model;

    model->set_owner(this);
    model->set_database(_database);
    if (this->get_type() == PhysicalObject::BULK)
      model->set_material(static_cast<const Material*>(this));
    model->set_simulator_id(simulator_id);
  }
}



void
PhysicalObject::init(void)
{
  if (!_is_initialized)
  {
    do_init();
    // free resources from database
    _database->close();
    _is_initialized = true;
  }
}



void
PhysicalObject::do_init(void)
{
  ModelMap::iterator it = models_begin();
  const ModelMap::const_iterator end = models_end();

  for ( ; it != end; ++it)
    (it->second)->init();
}



void
PhysicalObject::set_database(const Database& database)
{
  delete _database;
  _database = new Database(database);
}




