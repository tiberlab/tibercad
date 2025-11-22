// $Id$

#include "tibercad/physics/Alloy.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/io/Messages.h"
#include "tibercad/atomistic/CrystalDefs.h"


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

  // the fractions of possible sub-alloys
  double xa = -1.0;
  double xb = -1.0;

  if (get_options().find_option("comp_A"))
  {
    names.resize(2);
    names[0] = get_options()["comp_A"];
    names[1] = get_options()["comp_B"];
  }

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

  xa = get_database().get("x_A", xa);
  xb = get_database().get("x_B", xb);
  if (get_options().find_option("x_A"))
  {
    xa = get_options().get_option("x_A", xa);
    xb = get_options().get_option("x_B", xb);
  }

  get_database().set_alloy_mixing(mixing);

  {
    std::ostringstream os;
    os << get_name() << " is an alloy with components " <<
      names[0] << " and " << names[1] << " and molar fraction of " << names[0]
      << " "<< _molar_fraction << ".";
    Messages::info(os.str());
  }


  ModelOptions opts(get_options());
  opts.delete_option("comp_A");
  opts.delete_option("comp_B");
  opts.delete_option("x_A");
  opts.delete_option("x_B");
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
    modA->set_owner((it->second)->get_owner());
    //modA->set_database(new Database(_mat_A->get_database()));

    PhysicalModel* modB = static_cast<PhysicalModel*>((it->second)->copy());
    _mat_B->add_model(modB, it->first);
    modB->set_owner((it->second)->get_owner());
  }

  _mat_A->init();
  _mat_B->init();

  _cb_bloch_states = _mat_A->get_cb_bloch_functions();
  _vb_bloch_states = _mat_A->get_vb_bloch_functions();
  _cb_atomic_orbitals = _mat_A->get_cb_atomic_orbitals();
  _vb_atomic_orbitals = _mat_A->get_vb_atomic_orbitals();

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
  
  if (_specie_fraction.size() > 0) _specie_fraction.clear();
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
    _specie_fraction.resize((*it)->count_labels()+1);

    db.set_section("atomistic_structure");

    unsigned int ctr = 1;
    unsigned int n_species = db.get("n_basis_specie", 0);
    for (unsigned int i = 1; i <= n_species; i++)
    {
      std::stringstream out;
      out << i; 
      Specie tmp(db.get("specie_"+out.str(), "None"));
      _species.insert(tmp);

      std::string record = "n_" + out.str();
      unsigned int n = db.get(record.c_str(), 0);

      for (unsigned int j = 0; j < n; ++j, ++ctr)
      {
        double x = 0.0;
        if ( _specie_fraction[ctr].count(tmp))
        {
          x = _specie_fraction[ctr][tmp];
        }

        // x is used to sum up all fraction for the same specie.
        _specie_fraction[ctr][tmp] = x + mmap[*it];

        _crystal_type_map[ctr].insert(tmp);
      }
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


const Material*
Alloy::get_parent(std::pair<Specie, Specie> atom_pair) const
{
  const Material* mat = nullptr;
  const Material* mat1 = this->get_component_A();
  const Material* mat2 = this->get_component_B();
  if (mat1->has_specie(atom_pair.first) && mat1->has_specie(atom_pair.second))
  {
    // can take it from here
    mat = mat1;
  }
  else if (mat2->has_specie(atom_pair.first) && mat2->has_specie(atom_pair.second))
  {
    // can take it from here
    mat = mat2;
  }

  if ((mat != nullptr) && (mat->is_alloy()))
  {
    const Alloy* alloy = static_cast<const Alloy*>(mat);
    mat = alloy->get_parent(atom_pair);
  }

  return(mat);
}
