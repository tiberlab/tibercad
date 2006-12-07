// $Id$

#include "Material.h"
#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"

Database*
Material::_database;


Material::~Material(void)
{
  ModelMap::const_iterator it(_models.begin());
  const ModelMap::const_iterator end(_models.end());
  for ( ; it != end; ++it)
    delete it->second;

  _models.clear();

  delete _rotated_crystal;
}


void
Material::do_init(void)
{
  _database->set_material(_name);

  // first we set up RotatedCrystal because it will be
  // needed by others
  ModelOptions opts;
  if (get_options().find_option("a"))
  {
    opts["a"] = get_options()["a"];
    get_options().delete_option("a");
  }
  if (get_options().find_option("c"))
  {
    opts["c"] = get_options()["c"];
    get_options().delete_option("c");
  }
  if (get_options().find_option("x-growth-direction"))
  {
    opts["x-growth-direction"] = get_options()["x-growth-direction"];
    get_options().delete_option("x-growth-direction");
  }
  if (get_options().find_option("y-growth-direction"))
  {
    opts["y-growth-direction"] = get_options()["y-growth-direction"];
    get_options().delete_option("y-growth-direction");
  }
  if (get_options().find_option("z-growth-direction"))
  {
    opts["z-growth-direction"] = get_options()["z-growth-direction"];
    get_options().delete_option("z-growth-direction");
  }

  _rotated_crystal = RotatedCrystal::create(get_structure(), opts);
  _rotated_crystal->set_material(this);
  _rotated_crystal->init();

  ModelMap::iterator it = _models.begin();
  const ModelMap::const_iterator end = _models.end();

  for ( ; it != end; ++it)
  {
    (it->second)->set_material(this);

    (it->second)->init();
  }
}


void
Material::add_model(PhysicalModel* model, ID simulator_id)
{
  assert(model != NULL);
  assert(simulator_id != 0);
  
  ModelMap::iterator it = _models.find(simulator_id);
  if (it != _models.end())
  {
    delete it->second;
    it->second = model;
  }
  else
    _models[simulator_id] = model;

  model->set_simulator_id(simulator_id);
}


Material*
Material::create(const std::string& name)
{
  assert(_database != NULL);

  Material* mat = NULL;

  if (_database->is_alloy(name))
    mat = Alloy::create(name);
  else
    mat = new Material(name);

  std::cerr << "Created Material " << mat->get_name() << "\n";

  return mat;
}

Material*
Material::create(const std::string& name, const ModelOptions& options)
{
  Material* mat = create(name);

  if (mat != NULL)
  {
    mat->set_options(options);

    mat->_structure = mat->_options.get_option("structure", "zb");
    mat->_options.delete_option("structure");
  }

  return mat;
}
