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
  ModelOptions::const_submodel_iterator
      it(get_options().submodels_begin("Doping"));
  ModelOptions::const_submodel_iterator
      end(get_options().submodels_end("Doping"));
  for ( ; it != end; ++it)
  {
    const ModelOptions& opts = it->second;
    add_dopant(Dopant::create(opts.get_name(), opts));
  }
}




void
Material::do_init(void)
{
  setup_doping();

  PhysicalObject::do_init();
}



void
Material::do_preinit(void)
{
  ModelOptions opts;

  if (get_options().find_option("a"))
    opts["a"] = get_options()["a"];
  if (get_options().find_option("c"))
    opts["c"] = get_options()["c"];

  opts["x-growth-direction"] = get_options()["x-growth-direction"];
  opts["y-growth-direction"] = get_options()["y-growth-direction"];
  opts["z-growth-direction"] = get_options()["z-growth-direction"];

  // first we set up RotatedCrystal because it will be
  // needed by others
  _rotated_crystal = RotatedCrystal::create(get_structure(), opts);
  _rotated_crystal->set_owner(this);
  _rotated_crystal->init();
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
    mat->preinit();
  }

  return mat;
}



void
Material::preinit(void)
{

  // set the crystal structure at this point
  _structure = get_database().get("structure", "zb");

  bool hasx = get_options().find_option("x-growth-direction");
  bool hasy = get_options().find_option("y-growth-direction");
  bool hasz = get_options().find_option("z-growth-direction");

  // if one is given, all three should be provided
  if (((hasx || hasy) && !(hasx && hasy)) ||
      ((hasy || hasz) && !(hasy && hasz)))
    throw InitFailedException("You have to provide all growth directions "
        "in material " + get_name());

  bool use_defaults = !hasx;

  unsigned int dim = get_option("dimension", 4);

  // read or set default growth directions for wurtzite
  if (_structure == "wz")
  {
    if (dim == 4)
      Messages::warning("Material " + get_name() + " has no mesh dimension "
          "associated. Cannot guess crystal directions.");

    if (use_defaults)
    {
      switch (dim)
      {
        case 3:
          get_options()["x-growth-direction"] = "( 1,0,-1,0)";
          get_options()["y-growth-direction"] = "(-1,2,-1,0)";
          get_options()["z-growth-direction"] = "( 0,0, 0,1)";
          break;

        case 2:
          get_options()["z-growth-direction"] = "( 1,0,-1,0)";
          get_options()["x-growth-direction"] = "(-1,2,-1,0)";
          get_options()["y-growth-direction"] = "( 0,0, 0,1)";
          break;

        case 1:
        default:
          get_options()["y-growth-direction"] = "( 1,0,-1,0)";
          get_options()["z-growth-direction"] = "(-1,2,-1,0)";
          get_options()["x-growth-direction"] = "( 0,0, 0,1)";
      }
    }

  }
  else if (_structure == "zb")
  {
    if (use_defaults)
    {
      get_options()["x-growth-direction"] = "(1,0,0)";
      get_options()["y-growth-direction"] = "(0,1,0)";
      get_options()["z-growth-direction"] = "(0,0,1)";
    }
  }


  do_preinit();

  get_database().close();
}




void
Material::info(void) const
{
  std::ostringstream os;
  Messages m;
  m.indent();

  do_info();

  os << "x growth direction : " << get_option("x-growth-direction", "") << std::endl;
  os << "y growth direction : " << get_option("y-growth-direction", "") << std::endl;
  os << "z growth direction : " << get_option("z-growth-direction", "") << std::endl;

  m.info(os.str());
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




