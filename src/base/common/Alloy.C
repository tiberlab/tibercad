// $Id$

#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "PhysicalModel.h"
#include "Messages.h"


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
Alloy::do_preinit(void)
{
  std::string name_A, name_B;
  get_database().set_section("");
  get_database().get_alloy_components(name_A, name_B);

  Messages::debug(get_name() + " is an alloy with components " +
    name_A + " and " + name_B + ".");


  _mat_A = Material::create(name_A, get_options());
  _mat_B = Material::create(name_B, get_options());

  // to be sure we put the structure into the options
  _mat_A->set_structure(get_structure());
  _mat_B->set_structure(get_structure());
}
 

void
Alloy::do_init(void)
{
  // initialize the parent material
  //Material::do_init();

  setup_doping();

  _molar_fraction = get_options().get_option("x", 0.0);

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

  get_database().set_material(get_name());

  RotatedCrystal* crystal = static_cast<RotatedCrystal*>(
      _mat_A->get_rotated_crystal().copy());
  crystal->set_material(this);
  crystal->init_alloy(&_mat_A->get_rotated_crystal(),
      &_mat_B->get_rotated_crystal(), _molar_fraction);
  set_crystal(crystal);


  for (it = models_begin(); it != end; ++it)
  {
    (it->second)->init_alloy(_mat_A->get_model(it->first),
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


