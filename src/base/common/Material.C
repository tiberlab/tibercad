// $Id$

#include "Material.h"
#include "PhysicalModel.h"
#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Dopant.h"
#include "Messages.h"




Material::Material(const std::string& name,
    const ModelOptions& options, bool alloy)
  : PhysicalObject(BULK, options),
    _structure("zb"),
    _is_alloy(alloy),
    _rotated_crystal(NULL)
{
  set_name(name);
}



Material::~Material(void)
{
  clear_doping();

  destroy(_rotated_crystal);

  Messages::debug("Destroyed Material " + get_name());
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

  PhysicalObject::do_init();
}



void
Material::set_crystal(RotatedCrystal* crystal)
{
  destroy(_rotated_crystal);
  _rotated_crystal = crystal;
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
Material::create(const std::string& name, const ModelOptions& options)
{
  Material* mat = NULL;

  Database db(name, options.get_option("datafile", ""));
  db.set_section("");

  if (db.is_alloy())
    mat = Alloy::create(name, options);
  else
    mat = new Material(name, options);

  if (mat != NULL)
  {
    mat->set_database(db);

    Messages::debug("Created Material " + name +
      " (using parameter file " + db.get_data_file() + ")");

    mat->preinit();
  }

  return mat;
}



void
Material::preinit(void)
{

  // set the crystal structure at this point
  _structure = get_database().get("structure", "zb");

  do_preinit();

  get_database().close();
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




