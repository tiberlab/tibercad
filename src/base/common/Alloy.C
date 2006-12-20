// $Id$

#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"


inline
Alloy::Alloy(const std::string& name)
  : Material(name),
    _molar_fraction(0.0)
{
  get_database().get_alloy_components(name, _name_A, _name_B);
  std::cerr << name << " is an alloy with components " <<
    _name_A << " and " << _name_B << ".\n";
}


Alloy::~Alloy()
{
  ModelMap::const_iterator it(_models_A.begin());
  ModelMap::const_iterator end(_models_A.end());
  for ( ; it != end; ++it)
    PhysicalModelInterface::destroy(it->second);

  it = _models_B.begin();
  end = _models_B.end();
  for ( ; it != end; ++it)
    PhysicalModelInterface::destroy(it->second);

  _models_A.clear();
  _models_B.clear();

  PhysicalModelInterface::destroy(_cryst_A);
  PhysicalModelInterface::destroy(_cryst_B);
}


void
Alloy::do_init(void)
{
  // initialize the parent material
  Material::do_init();

  _molar_fraction = get_options().get_option("x", 0.0);

  // copy and initialize the models of the components
  ModelMap::iterator it(models_begin());
  ModelMap::const_iterator end(models_end());

  for ( ; it != end; ++it)
  {
    PhysicalModel* modA = static_cast<PhysicalModel*>((it->second)->copy());
    modA->set_material(this);
    get_database().set_material(_name_A);
    modA->init();
    _models_A[it->first] = modA;
    
    PhysicalModel* modB = static_cast<PhysicalModel*>((it->second)->copy());
    modB->set_material(this);
    get_database().set_material(_name_B);
    modB->init();
    _models_B[it->first] = modB;

    // build VCA
    get_database().set_material(get_name());
    (it->second)->build_alloy(modA, modB, _molar_fraction);
  }

  _cryst_A = static_cast<RotatedCrystal*>(get_rotated_crystal().copy());
  _cryst_A->set_material(this);
  get_database().set_material(_name_A);
  _cryst_A->init();

  _cryst_B = static_cast<RotatedCrystal*>(get_rotated_crystal().copy());
  _cryst_B->set_material(this);
  get_database().set_material(_name_B);
  _cryst_B->init();

  get_database().set_material(get_name());
  get_crystal().build_alloy(_cryst_A, _cryst_B, _molar_fraction);
}


Alloy*
Alloy::create(const std::string& name)
{
  Alloy* mat = NULL;

  mat = new Alloy(name);

  return mat;
}


