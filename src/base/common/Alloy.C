// $Id$

#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"


inline
Alloy::Alloy(const std::string& name)
  : Material(name),
    _molar_fraction(0.0),
    _mat_A(NULL),
    _mat_B(NULL)
{
  _is_alloy = true;
}


Alloy::~Alloy()
{
  delete _mat_A;
  delete _mat_B;
}


void
Alloy::do_init(void)
{
  // initialize the parent material
  Material::do_init();

  _molar_fraction = get_options().get_option("x", 0.0);

  std::string name_A, name_B;
  get_database().get_alloy_components(get_name(), name_A, name_B);
#ifdef DEBUG
  std::cout << get_name() << " is an alloy with components " <<
    name_A << " and " << name_B << ".\n";
#endif
  
  _mat_A = Material::create(name_A, get_options());
  _mat_B = Material::create(name_B, get_options());
  
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

  get_database().set_material(get_name());

  get_crystal().build_alloy(&_mat_A->get_rotated_crystal(),
      &_mat_B->get_rotated_crystal(), _molar_fraction);

  for (it = models_begin(); it != end; ++it)
  {
    (it->second)->build_alloy(_mat_A->get_model(it->first),
                              _mat_B->get_model(it->first), _molar_fraction);
  }
}


Alloy*
Alloy::create(const std::string& name)
{
  Alloy* mat = NULL;

  mat = new Alloy(name);

  return mat;
}


