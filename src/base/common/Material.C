// $Id$

#include "Material.h"
#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Dopant.h"

#include "getpot.h"

Database*
Material::_database;



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
Material::do_init(void)
{
  _database->set_material(_name);

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

  // first we set up RotatedCrystal because it will be
  // needed by others
  _rotated_crystal = RotatedCrystal::create(get_structure(), opts);
  _rotated_crystal->set_material(this);
  _rotated_crystal->init();


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
    PhysicalModelInterface::destroy(it->second);
    it->second = model;
  }
  else
    _models[simulator_id] = model;

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
  std::cout << "Created Material " << mat->get_name() << 
    " (using parameter file " << _database->get_data_file() << ")" <<
    std::endl;

  return mat;
}




Material*
Material::create(const std::string& name, const ModelOptions& options)
{

  Material* mat = create(name);

  if (mat != NULL)
  {
    mat->set_options(options);

    _database->set_material(mat->get_name());
    GetPot data(_database->get_data_file());
    mat->_structure = data("structure", "zb");

    // set the crystal structure at this point
    mat->_structure = mat->_options.get_option("structure", mat->_structure);
    mat->_options.delete_option("structure");
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



