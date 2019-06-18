// $Id$

#include "Material.h"
#include "PhysicalModel.h"
#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "Dopant.h"
#include "Messages.h"
#include "CrystalDefs.h"



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
  if (get_options().find_option("b"))
    opts["b"] = get_options()["b"];
  if (get_options().find_option("c"))
    opts["c"] = get_options()["c"];

  if (get_options().find_option("alpha"))
    opts["alpha"] = get_options()["alpha"];
  if (get_options().find_option("beta"))
    opts["beta"] = get_options()["beta"];
  if (get_options().find_option("gamma"))
    opts["gamma"] = get_options()["gamma"];

  opts["x-growth-direction"] = get_options()["x-growth-direction"];
  opts["y-growth-direction"] = get_options()["y-growth-direction"];
  opts["z-growth-direction"] = get_options()["z-growth-direction"];

  // first we set up RotatedCrystal because it will be
  // needed by others
  _rotated_crystal = RotatedCrystal::create(this, opts);
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

  // first check if an alloy is specified
  if (options.has_submodel("alloy"))
  {
    ModelOptions alloy_opts(options.submodels_begin("alloy")->second);

    ModelOptions::submodel_iterator comp_it(alloy_opts.submodels_begin("component"));
    const ModelOptions::submodel_iterator comp_end(alloy_opts.submodels_end("component"));

    if (std::distance(comp_it, comp_end) != 2)
      throw InitFailedException("Currently alloys with more than two components "
          "cannot be defined.");

    //for ( ; comp_it != comp_end; ++comp_it)
    //{
    //
    //}

    std::string comp_A(comp_it->second.get_name());
    double x_A = comp_it->second.get_option("x", 0.5);
    ++comp_it;
    std::string comp_B(comp_it->second.get_name());
    double x_B = comp_it->second.get_option("x", 0.5);

    // we use this to disable database files
    alloy_opts["datafile"] = alloy_opts.get_option("datafile", "");

    mat = Alloy::create(alloy_opts.get_name(), alloy_opts);

  }
  else
  {
    Database db(name, options.get_option("datafile", ""));
    db.set_section("");

    if (db.is_alloy())
      mat = Alloy::create(name, options);
    else
      mat = new Material(name, options);

  }

  if (mat != NULL)
  {
    mat->preinit();
  }

  return mat;
}



void
Material::preinit(void)
{
  // set the database
  Database db(get_name(), get_options().get_option("datafile", ""));
  set_database(db);

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
  
  if (dim == 4)
      Messages::warning("Material " + get_name() + " has no mesh dimension "
          "associated. Cannot guess crystal directions.");

  if (use_defaults)
  {
    get_options()["x-growth-direction"] = "(1,0,0)";
    get_options()["y-growth-direction"] = "(0,1,0)";
    get_options()["z-growth-direction"] = "(0,0,1)";
  }

  // read or set default growth directions for wurtzite
  if (_structure == "wz")
  {
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
          break;
      }
    }

  }


  do_preinit();

  fill_species();

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



void
Material::fill_species(void)
{
  //Read and store species belonging to the material.
  //Useful for assigning components to atoms in random alloys
  get_database().set_section("atomistic_structure");
  unsigned int n_species = get_database().get("n_basis_specie", 0);
  for (unsigned int i = 1; i <= n_species; i++)
  {
    std::stringstream out;
    out << i; 
    Specie tmp(get_database().get("specie_" + out.str(), "None"));
    _species.insert(tmp);
    _crystal_type_map[i].insert(tmp);
  }
}


bool
Material::has_specie(Specie sp) const
{
  bool has_it = false;

  if (this->is_alloy())
  {
    const Alloy* alloy = static_cast<const Alloy*>(this);
    has_it = alloy->get_component_A()->has_specie(sp);

    if (!has_it)
      has_it = alloy->get_component_B()->has_specie(sp);
  }
  else
    has_it = (_species.find(sp) != _species.end());

  return has_it;
}


bool
Material::is_specie(Specie sp, unsigned int label) const
{
  if (!has_specie(sp))
  {
    Messages::error("Error in is_specie(): specie " +
                    sp.get_string() + "not found in material " + get_name());
  }
  else
  {
    return ( ((_crystal_type_map.find(label))->second).count(sp) > 0 ); 
  }

}

unsigned int
Material::get_label(Specie sp) const
{
  if (!has_specie(sp))
  {
    // TODO no magic numbers, please...
    return 255;
  }
  else
  {
    std::map<unsigned int, std::set<Specie>>::const_iterator
      it = _crystal_type_map.begin();
    
    for ( ; it != _crystal_type_map.end(); ++it)
      if ( (it->second).count(sp) > 0 ) return it->first; 
    
  }
}

void
Material::print_species(void) const
{
  std::cout<<"Species: ";
  std::set<Specie>::iterator sp = _species.begin();
  for( ; sp != _species.end(); ++sp)
    std::cout<< (*sp)<<"  ";
  std::cout<<std::endl;
  
  for (unsigned int i=1; i<= _crystal_type_map.size(); i++)
  {
    std::set<Specie>::iterator it =  (_crystal_type_map.find(i)->second).begin();
    std::cout<<i;
    
    for( ; it != (_crystal_type_map.find(i)->second).end(); ++it)
    {
      std::cout<<"  "<< (*it);
    }
    std::cout<<std::endl;
  }
}


