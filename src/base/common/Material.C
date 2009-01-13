// $Id$

#include "Material.h"
#include "PhysicalModel.h"
#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Dopant.h"

#include "getpot.h"

Database*
Material::_database;


Material::Material(const std::string& name)
  : _is_alloy(false),
    _name(name),
    _structure("zb"), 
    _rotated_crystal(NULL),
    _is_initialized(false)
{
}



Material::~Material(void)
{
  ModelMap::const_iterator it(_models.begin());
  const ModelMap::const_iterator end(_models.end());
  for ( ; it != end; ++it)
    PhysicalModelInterface::destroy(it->second);

  _models.clear();

  clear_doping();

  PhysicalModelInterface::destroy(_rotated_crystal);
}



void
Material::setup_doping(void)
{
  // now the doping
  double doping = get_options().get_option("doping", 0.0);
  if (doping > 0.0)
  {
    double level = get_options().get_option("doping_level", 0.025);
    int g = get_options().get_option("doping_degen", 2);
    // allow simplified name
    g = get_options().get_option("g", g);
    Dopant::DopingType type = Dopant::N_TYPE;
    const std::string& doptype = get_options().get_option("doping_type", "");
    if (doptype == "acceptor")
      type = Dopant::P_TYPE;
    
    add_dopant(new Dopant(doping, level, g, type));
  }
}




void
Material::do_init(void)
{
  ModelOptions opts;

  if (get_options().find_option("a"))
    opts["a"] = get_options()["a"];
  if (get_options().find_option("c"))
    opts["c"] = get_options()["c"];
  if (get_options().find_option("x-growth-direction"))
    opts["x-growth-direction"] = get_options()["x-growth-direction"];
  if (get_options().find_option("y-growth-direction"))
    opts["y-growth-direction"] = get_options()["y-growth-direction"];
  if (get_options().find_option("z-growth-direction"))
    opts["z-growth-direction"] = get_options()["z-growth-direction"];


  // first we set up RotatedCrystal because it will be
  // needed by others
  _rotated_crystal = RotatedCrystal::create(get_structure(), opts);
  _rotated_crystal->set_material(this);
  _rotated_crystal->init();

  setup_doping();

  ModelMap::iterator it = _models.begin();
  const ModelMap::const_iterator end = _models.end();

  for ( ; it != end; ++it)
    (it->second)->init();
}


void
Material::set_crystal(RotatedCrystal* crystal)
{
  PhysicalModelInterface::destroy(_rotated_crystal);
  _rotated_crystal = crystal;
}



void
Material::add_model(PhysicalModel* model, ID simulator_id)
{
  assert(model != NULL);
  assert(simulator_id != 0);
  
  ModelMap::iterator it = _models.find(simulator_id);
  if (it != _models.end())
  {
    PhysicalModelInterface::destroy(it->second);
    it->second = model;
  }
  else
    _models[simulator_id] = model;

  model->set_material(this);
  model->set_simulator_id(simulator_id);
}



void
Material::add_dopant(Dopant* dopant)
{
  if (dopant != NULL)
  {
    if (dopant->get_type() == Dopant::N_TYPE)
      _donors.insert(dopant);
    else
      _acceptors.insert(dopant);
  }
    
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

  _database->set_material(name);

  return mat;
}



Material*
Material::create(const std::string& name, const ModelOptions& options)
{
  Material* mat = NULL;

  _database->set_material(name, options.get_option("datafile", ""));
  _database->set_section("");

  if (_database->is_alloy(name))
    mat = Alloy::create(name);
  else
    mat = new Material(name);

  if (mat != NULL)
  {
    mat->set_options(options);

    // set the crystal structure at this point
    mat->_structure = mat->get_database().get("structure", "zb");

    std::cout << "Created Material " << mat->get_name() << 
      " (using parameter file " << _database->get_data_file() << ")" <<
      std::endl;

    mat->preinit();
  }

  return mat;
}



double
Material::get_total_donor_density(void) const
{
  double Nd = 0;
  dopant_iterator it = _donors.begin();
  dopant_iterator end = _donors.end();
  for ( ; it != end; ++it)
    Nd += (*it)->get_doping_density();

  return Nd;
}



double
Material::get_total_acceptor_density(void) const
{
  double Na = 0;
  dopant_iterator it = _acceptors.begin();
  dopant_iterator end = _acceptors.end();
  for ( ; it != end; ++it)
    Na += (*it)->get_doping_density();

  return Na;
}



void
Material::clear_doping(void)
{
  dopant_iterator it = _donors.begin();
  dopant_iterator end = _donors.end();
  for ( ; it != end; ++it)
    delete (*it);

  it = _acceptors.begin();
  end = _acceptors.end();
  for ( ; it != end; ++it)
    delete (*it);

  _donors.clear();
  _acceptors.clear();
}



const Database&
Material::get_database(void) const
{
  assert(_database != NULL);
  _database->set_material(get_name(),
      get_options().get_option("datafile", ""));
  return *_database;
}


Database&
Material::get_database(void)
{
  assert(_database != NULL);
  _database->set_material(get_name(),
      get_options().get_option("datafile", ""));
  return *_database;
}


