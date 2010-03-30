// $Id$

#include "Alloy.h"
#include "Database.h"
#include "RotatedCrystal.h"
#include "PhysicalModel.h"
#include "Messages.h"


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

#ifdef DEBUG
  std::ostringstream os;
  os << get_name() << " is an alloy with components " <<
    names[0] << " and " << names[1] << " and molar fraction of " << names[0]
    << " is " << _molar_fraction << ".";
  Messages::debug(os.str());
#endif



  _mat_A = Material::create(names[0], get_options());
  _mat_B = Material::create(names[1], get_options());


  // TODO for now alloy has two components
  assert(get_database().get_number_of_components() == 2);

  std::vector<double> fracs(2);
  fracs[0] = _molar_fraction;
  fracs[1] = 1.0 - _molar_fraction;
  get_database().set_alloy_composition(fracs);

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

    PhysicalModel* modB = static_cast<PhysicalModel*>((it->second)->copy());
    _mat_B->add_model(modB, it->first);
  }

  _mat_A->init();
  _mat_B->init();

  //
  // build VCA of the models
  //

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
Alloy::create(const std::string& name, const ModelOptions& options)
{
  Alloy* mat = NULL;

  mat = new Alloy(name, options);

  return mat;
}


