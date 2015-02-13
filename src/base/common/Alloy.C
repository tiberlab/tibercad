// $Id$

#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "PhysicalModel.h"
#include "Messages.h"
#include "CrystalDefs.h"


Alloy::Alloy(const std::string& name, const ModelOptions& options)
  : Material(name, options, true),
    _molar_fraction(0.0),
    _mat_A(NULL),
    _mat_B(NULL)
{
}


Alloy::~Alloy()
{
  delete _mat_A;
  delete _mat_B;
}



void
Alloy::do_preinit(void)
{
  std::vector<std::string> names;
  get_database().set_section("");
  get_database().get_components(names);

  // for now, alloys can have only two components!
  assert(names.size() == 2);

  _molar_fraction = get_options().get_option("x", 0.0);

  std::vector<double> fracs(2);
  fracs[0] = _molar_fraction;
  fracs[1] = 1.0 - _molar_fraction;
  get_database().set_alloy_composition(fracs);

  // the components may be alloys
  Database::AlloyMixing mixing = get_database().get_alloy_mixing();
  get_database().set_alloy_mixing(Database::NONE);
  double xa = get_database().get("x_A", -1.0);
  double xb = get_database().get("x_B", -1.0);
  get_database().set_alloy_mixing(mixing);

  {
    std::ostringstream os;
    os << get_name() << " is an alloy with components " <<
      names[0] << " and " << names[1] << " and molar fraction of " << names[0]
      << " "<< _molar_fraction << ".";
    Messages::info(os.str());
  }


  ModelOptions opts(get_options());
  if (xa >= 0)
    opts.set_option("x", xa);
  _mat_A = Material::create(names[0], opts);

  opts.delete_option("x");
  if (xb >= 0)
    opts.set_option("x", xb);
  _mat_B = Material::create(names[1], opts);

  // if components are alloys, we have to set correctly their alloy fractions
  if (xa >= 0)
  {
    fracs[0] = xa;
    fracs[1] = 1.0 - xa;
    get_database().get_component_database(0).set_alloy_composition(fracs);
  }
  if (xb >= 0)
  {
    fracs[0] = xb;
    fracs[1] = 1.0 - xb;
    get_database().get_component_database(1).set_alloy_composition(fracs);
  }

  _mat_A->set_database(get_database().get_component_database(0));
  _mat_B->set_database(get_database().get_component_database(1));

  // a sanity check on the crystal structure
  if ((_mat_A->get_structure() != get_structure()) ||
      (_mat_B->get_structure() != get_structure()))
  {
    std::ostringstream os;
    os << "The crystal structures of the Alloy " << get_name()
       << " and its components are inconsistent.";
    throw InitFailedException(os.str());
  }

  RotatedCrystal* crystal = static_cast<RotatedCrystal*>(
      _mat_A->get_rotated_crystal().copy());
  crystal->set_owner(this);
  crystal->init_alloy(&_mat_A->get_rotated_crystal(),
      &_mat_B->get_rotated_crystal(), _molar_fraction);
  set_crystal(crystal);



}


void
Alloy::do_init(void)
{

  setup_doping();

  assert((_mat_A != NULL) && (_mat_B != NULL));

  // copy and initialize the models of the components
  ModelMap::iterator it(models_begin());
  ModelMap::const_iterator end(models_end());

  for ( ; it != end; ++it)
  {
    PhysicalModel* modA = static_cast<PhysicalModel*>((it->second)->copy());
    _mat_A->add_model(modA, it->first);

    PhysicalModel* modB = static_cast<PhysicalModel*>((it->second)->copy());
    _mat_B->add_model(modB, it->first);
  }

  _mat_A->init();
  _mat_B->init();

  //
  // build VCA of the models
  //


  for (it = models_begin(); it != end; ++it)
  {
    (it->second)->init_alloy(_mat_A->get_model(it->first),
                             _mat_B->get_model(it->first), _molar_fraction);
  }
}


void
Alloy::do_info(void) const
{
  std::ostringstream os;
  os << "alloy with " << _molar_fraction << " " << _mat_A->get_name()
      << ", " << (1 - _molar_fraction) << " " << _mat_B->get_name();
  Messages::info(os.str());
}

void
Alloy::fill_species(void) 
{
  std::set<const Material*> materials;
  std::map<const Material*, double> mmap;
  
  if (_specie_fraction.size()>0) _specie_fraction.clear();
  _specie_fraction.resize(3); 
  _species.clear();

  // Ga(x)In(1-x)N =>  GaN(x) InN(1-x)
  if (_mat_A->is_alloy()) // this is quaternary case
  {
    const Alloy* alloy = dynamic_cast<const Alloy*>(_mat_A);
    materials.insert(alloy->get_component_A());
    materials.insert(alloy->get_component_B());    
    mmap[alloy->get_component_A()] = 
      (alloy->get_molar_fraction()) * _molar_fraction;
    mmap[alloy->get_component_B()] = 
      (1.0-alloy->get_molar_fraction()) * _molar_fraction;
  }
  else
  {
    materials.insert(_mat_A);
    mmap[_mat_A] = _molar_fraction;
  } 

  if (_mat_B->is_alloy()) // this is quaternary case
  {
    const Alloy* alloy = dynamic_cast<const Alloy*>(_mat_B);
    materials.insert(alloy->get_component_A());
    materials.insert(alloy->get_component_B()); 
    mmap[alloy->get_component_A()] = 
      (alloy->get_molar_fraction()) * (1.0-_molar_fraction);
    mmap[alloy->get_component_B()] = 
      (1.0-alloy->get_molar_fraction()) * (1.0-_molar_fraction);
  }
  else
  {
    materials.insert(_mat_B);
    mmap[_mat_B] = 1.0 - _molar_fraction;
  }
   
  std::set<const Material*>::const_iterator it = materials.begin();
  for( ; it != materials.end(); ++it)
  {
    const Database& db = (*it)->get_database();
    db.set_section("atomistic_structure");
    unsigned int n_species = db.get("n_basis_specie", 0);
    for (unsigned int i = 1; i <= n_species; i++)
    {
      std::stringstream out;
      out << i; 
      Specie tmp(db.get("specie_"+out.str(), "None"));
      _species.insert(tmp);

      if (_specie_fraction.size() < i) _specie_fraction.resize(i+1); 

      double x=0.0;
      if ( _specie_fraction[i].count(tmp)){ x=_specie_fraction[i][tmp]; }

      // x is used to sum up all fraction for the same specie.
      _specie_fraction[i][tmp] = x+mmap[*it];

      _crystal_type_map[i].insert(tmp);
    }
 

  }
  
  /*
  std::map<Specie,double>::iterator spit = _specie_fraction[1].begin();
  for( ; spit != _specie_fraction[1].end(); ++spit)
  {    
    std::cout<<"(Alloy) "<<spit->first<<" "<<spit->second<<std::endl;
  }
  spit = _specie_fraction[2].begin();
  for( ; spit != _specie_fraction[2].end(); ++spit)
  {    
    std::cout<<"(Alloy) "<<spit->first<<" "<<spit->second<<std::endl;
  }
  */
  
}

bool
Alloy::is_mutable(unsigned int i) const
{
  return _specie_fraction[i].size() > 1;
}

Alloy*
Alloy::create(const std::string& name, const ModelOptions& options)
{
  Alloy* mat = NULL;

  mat = new Alloy(name, options);

  return mat;
}

/*
bool
Alloy::is_anion(Specie sp) const
{
  if (!((_mat_A->has_specie(sp)) || (_mat_B->has_specie(sp))))
    {
      Messages::error("Error in is_anion: specie not defined for any parent material");
    }
    else
    {
      return 
      (CrystalDefs::is_anion(_mat_A->get_name(), sp) || 
       CrystalDefs::is_anion(_mat_B->get_name(), sp));
    }
}

bool
Alloy::is_cation(Specie sp) const
{
  if (!((_mat_A->has_specie(sp)) || (_mat_B->has_specie(sp))))
    {
      Messages::error("Error in is_cation: specie" + sp.get_string() + " not defined for any parent material");
    }
    else
    {
      return 
      (CrystalDefs::is_cation(_mat_A->get_name(), sp) || 
       CrystalDefs::is_cation(_mat_B->get_name(), sp));
    }
}
*/
