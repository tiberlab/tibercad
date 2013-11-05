// $Id$

#include "tiber_config.h"
#include "PhysicalModelInterface.h"
#include "MaterialBoundary.h"
#include "Material.h"
#include "Variable.h"
#include "Database.h"

#include "Utils.h"
#include "Trap.h"
#include "ParticleDensity.h"

#include "ZbStiffness.h"
#include "WzStiffness.h"
#include "ZbPiezoelectricity.h"
#include "WzPiezoelectricity.h"
#include "ZbRotatedCrystal.h"
#include "WzRotatedCrystal.h"
#include "MacrostrainModel.h"
#include "ZbSemiconductor.h"
#include "WzSemiconductor.h"
#include "ZbDDsemiconductor.h"
#include "WzDDsemiconductor.h"

#include "SBZbCondBandBulkHamiltonian.h"
#include "SBWzCondBandBulkHamiltonian.h"
#include "SBuserHamiltonian.h"
#include "KPbulkHamiltonian.h"
#include "EFAbulkModel.h"

//#include  "DftbModel.h"
//#include  "EtbModel.h"


//#include  "ZbOptDielectricConstant.h"
//#include  "WzOptDielectricConstant.h"

//#include "PhononModel.h"
//#include "ZbFreeDynamicalMatrix.h"
//#include "ZbStrainDynamicalMatrix.h"
//#include "ZbRamanTensor.h"



#include "Messages.h"

#include <string>



using namespace std;


map<const string, ID>
PhysicalModelInterface::_model_ids;



PhysicalModelInterface::~PhysicalModelInterface(void)
{
  SubmodelIterator it(submodels_begin());
  const SubmodelIterator end(submodels_end());
  for ( ; it != end; ++it)
    destroy(it->second);
}



PhysicalModelInterface*
PhysicalModelInterface::_create(const string& name,
    const PhysicalObject* owner,
    const ModelOptions& options,
    const string& module)
{

  PhysicalModelInterface* mod = NULL;

  if (name == "stiffness_zb")
    mod = ZbStiffness::create(options);
  else if (name == "stiffness_wz")
    mod = WzStiffness::create(options);
  else if (name == "piezo_zb")
    mod = ZbPiezoelectricity::create(options);
  else if (name == "piezo_wz")
    mod = WzPiezoelectricity::create(options);
  else if (name == "cryst_zb")
    mod = ZbRotatedCrystal::create(options);
  else if (name == "cryst_wz")
    mod = WzRotatedCrystal::create(options);
  else if (name == "macrostrain")
    mod = MacrostrainModel::create(options);
  else if (name == "semicond_zb")
    mod = ZbSemiconductor::create(options);
  else if (name == "semicond_wz")
    mod = WzSemiconductor::create(options);
  else if (name == "quantum_kp")
    mod = KPbulkHamiltonian::create(options);
  else if (name == "quantum_cond_band_zb")
    mod = SBZbCondBandBulkHamiltonian::create(options);
  else if (name == "quantum_cond_band_wz")
    mod = SBWzCondBandBulkHamiltonian::create(options);
  else if (name == "quantum_user")
    mod = SBuserHamiltonian::create(options);
  else if (name == "DDsemicond_zb")
    mod = ZbDDsemiconductor::create(options);
  else if (name == "DDsemicond_wz")
    mod = WzDDsemiconductor::create(options);
  else if (name == "EFAmodel")
    mod = EFAbulkModel::create(options);
  //else if (name == "dftb")
  //  mod = DftbModel::create(options);
  //else if (name == "etb")
  //  mod = ETBModel::create(options);
  //else if (name == "opt_dielectric_constant_zb")
  //  mod = ZbOptDielectricConstant::create(options);
  //else if (name == "opt_dielectric_constant_wz")
  //  mod = WzOptDielectricConstant::create(options);
  //else if (name == "maxwell")
  //  mod = MaxwellPhysicalModel::create(options);
  //else if  (name == "phonon")
  //  mod = PhononModel::create(options);
  //else if  (name == "free_dynamical_matrix_zb")
  //  mod = ZbFreeDynamicalMatrix::create(options);
  //else if  (name == "strain_dependent_zb")
  //  mod = ZbStrainDynamicalMatrix::create(options);
  //else if  (name == "raman_tensor_zb")
  //  mod = ZbRamanTensor::create(options);
  //else if  (name == "piezoelectric_model_zb")
  //  mod = ZbPiezoelectricModel::create(options);
  //else if  (name == "piezoelectric_model_wz")
  //  mod = WzPiezoelectricModel::create(options);
  else if (name == "trap")
    mod = Trap::create(options);
  else if (name == "particle_density")
    mod = ParticleDensity::create(options);


  if (mod == NULL)
  {
    // first try without module directory
    if ((mod = create_from_library<PhysicalModelInterface>(name, options, owner)) == 0)
    {
      if (module.size() >  0)
      {
        mod = create_from_library<PhysicalModelInterface>(
            module + "/" + name, options, owner);
      }
    }
  }

  if (mod != NULL)
  {
    _register_model(mod);

    // we let it know what's its identifier
    mod->_set_type(name);

    mod->_set_module_name(module);
    mod->set_owner(owner);

    //! set the name
    // 2007-08-17
    //    we don't set anymore a default name
    //string defaultname = mod->get_default_name();
    //string defaultname = "";
    //mod->set_name(mod->get_options().get_option("name", defaultname));
    //mod->get_options().delete_option("name");

    ostringstream os;
    os << "Added model (ID = " << mod->get_id() <<
      " name = " << mod->get_name() << " type = " <<
       mod->get_type() << " module = " << mod->get_module_name() << ")";
    Messages::debug(os.str());
  }

  return mod;
}



PhysicalModelInterface*
PhysicalModelInterface::_create(create_t create_fnc, destroy_t destroy_fnc,
    const PhysicalObject* owner,
    const ModelOptions& options,
    const string& module)
{
  PhysicalModelInterface* mod = dynamic_cast<PhysicalModelInterface*>(
      create_from_function(create_fnc, destroy_fnc, options, owner));

  if (mod != NULL)
  {
    _register_model(mod);

    // we let it know what's its identifier
    mod->_set_type(options.get_option("type", ""));

    mod->_set_module_name(module);
    mod->set_owner(owner);

    //! set the name
    // 2007-08-17
    //    we don't set anymore a default name
    //string defaultname = mod->get_default_name();
    //string defaultname = "";
    //mod->set_name(mod->get_options().get_option("name", defaultname));
    mod->get_options().delete_option("name");
#ifdef DEBUG
    cerr << "Add model (ID = " << mod->get_id() <<
      " name = " << mod->get_name() << " type_id = " <<
       mod->get_type() << ")\n";
#endif
  }

  return mod;
}






void
PhysicalModelInterface::_register_model(
    PhysicalModelInterface* model)
{
  const string name = typeid(*model).name();
  model_id_iterator it = _model_ids.find(name);

  ID id;

  // set the model ID (create a new one if the model didn't exist yet)
  if (it == _model_ids.end())
  {
    id = _model_ids.size() + 1;
    _model_ids[name] = id;
  }
  else
    id = it->second;

  model->_id = id;
}



void
PhysicalModelInterface::set_owner(const PhysicalObject* owner)
{
  _owner = owner;
  if ((_owner != NULL) && (_owner->get_type() == PhysicalObject::BULK))
    _bulk_material = static_cast<const Material*>(_owner);
  else
    _bulk_material = NULL;
}



PhysicalModelInterface*
PhysicalModelInterface::copy(void) const
{
  PhysicalModelInterface* new_copy = NULL;

  // this is safe
  new_copy = this->create_new();

  if (new_copy != NULL)
  {
    new_copy->_id = _id;
    new_copy->set_owner(_owner);
    new_copy->set_material(_bulk_material);
    new_copy->set_name(get_name());
    new_copy->_set_type(get_type());
    new_copy->_simulator_id = _simulator_id;
    new_copy->_set_module_name(get_module_name());

    new_copy->copy_from(this);
  }

  return new_copy;
}



PhysicalModelInterface*
PhysicalModelInterface::create_new(void) const
{
  PhysicalModelInterface* pmi = create_from_object(this, get_owner());
  if (pmi == NULL)
  {
    ostringstream os;
    //os << "Model " << get_name() << " cannot create a new instance of "
    os << "Model cannot create a new instance of "
        "the same type as the method \"create_new()\" is not "
        "reimplemented.";
    throw ModelErrorException(os.str());
  }

  return pmi;
}



string
PhysicalModelInterface::get_default_name(void) const
{
  return Utils::extract_typename(typeid(*this));
}




const Database&
PhysicalModelInterface::get_database(void)
{
  return _owner->get_database();
}



void
PhysicalModelInterface::override_parameter_string(const std::string& name,
        std::string& s) const
{
  const PhysicalObject* mat = get_owner();

  if (mat != NULL)
  {
    std::string code;

    SimulationInterface* sim =
      SimulationInterface::get_simulation(get_simulator_id());

    // first ask for the plain name
    s = mat->get_options().get_option(name, s);

    // first override
    if (sim != NULL)
    {
      code = sim->get_name() + ".";
      s = mat->get_options().get_option(code + name, s);
    }

    // second override
    if (get_name() != "")
    {
      code += get_name() + ".";
      s = mat->get_options().get_option(code + name, s);
    }
  }
}






void
PhysicalModelInterface::init(void)
{
  //
  // we may arrive here, although the model is associated to an alloy
  // in that case, we defer to the init_interface() routine
  //
  if (get_owner()->get_type() == PhysicalObject::BOUNDARY)
  {
    const MaterialBoundary* mb = static_cast<const MaterialBoundary*>(get_owner());
    init_interface(mb->get_material_A(), mb->get_material_B());

    // NOTE we immediately return here
    return;
  }

  read_database();

  _create_submodels();
  Messages::debug("init() of " + get_material()->get_name());

  // Now we initialize all "official" submodels
  SubmodelIterator smit(submodels_begin());
  const SubmodelIterator smend(submodels_end());
  for ( ; smit != smend; ++smit)
  {
    PhysicalModelInterface* pm = smit->second;
    //pm->_simulator_id = _simulator_id;
    Messages::debug("Initializing " + smit->first + " ...");
    pm->init();
  }

  do_init();

  // dummy read
  get_option("crystal_structure", "");
  get_option("regions", "");
  get_option("type", "");
  get_option("name", "");

  get_options().check_unused(1);
}




void
PhysicalModelInterface::reinit(void)
{
  // reinit submodels
  SubmodelIterator smit(submodels_begin());
  const SubmodelIterator smend(submodels_end());
  for ( ; smit != smend; ++smit)
    smit->second->reinit();

  do_reinit();
}





void
PhysicalModelInterface::init_interface(const Material* comp_A,
    const Material* comp_B)
{
  // sometimes this is not a good idea, and it even may not make
  // sense at all in general
  //read_database();

  read_interface_database();

  // setup the submodels
  _create_submodels();
  Messages::debug("init() of " + get_owner()->get_name());

  SubmodelIterator it(submodels_begin());
  const SubmodelIterator end(submodels_end());

  for ( ; it != end; ++it)
    it->second->init_interface(comp_A, comp_B);

  // some models might treat interfaces in a special way
  do_init_interface(comp_A, comp_B);

  // dummy read
  get_option("regions", "");
  get_option("type", "");
  get_option("name", "");

  get_options().check_unused(1);
}


void
PhysicalModelInterface::init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  assert(typeid(*comp_A) == typeid(*comp_B));

  //
  // NOTE:
  //   The current approach relies on the strong assumption that all models
  //   in all alloy components are ordered exactly the same way. This had
  //   better be changed in the future!

  read_database();

  // some models might treat alloys in a special way
  // disable alloy mixing
  // This const cast is very ugly, better would be to not touch the database at all
  Database::AlloyMixing mixing = get_database().get_alloy_mixing();
  const_cast<PhysicalObject*>(_owner)->get_database().set_alloy_mixing(Database::NONE);
  read_database_alloy();
  const_cast<PhysicalObject*>(_owner)->get_database().set_alloy_mixing(mixing);


  // setup the submodels
  _create_submodels();

  SubmodelIterator it(submodels_begin());
  ConstSubmodelIterator it_A(comp_A->submodels_begin());
  ConstSubmodelIterator it_B(comp_B->submodels_begin());
  const SubmodelIterator end(submodels_end());
  for ( ; it != end; ++it, ++it_A, ++it_B)
  {
    PhysicalModelInterface* pm = it->second;
    pm->init_alloy(it_A->second, it_B->second, xa);
  }

  // some models might treat alloys in a special way
  do_init_alloy(comp_A, comp_B, xa);

  do_init();

  // dummy read
  get_option("crystal_structure", "");
  get_option("regions", "");
  get_option("type", "");
  get_option("name", "");

  get_options().check_unused(1);
}


void
PhysicalModelInterface::set_material(const Material* mat)
{
  _bulk_material = mat;

  SubmodelIterator it(_submodels.begin());
  const SubmodelIterator end(_submodels.end());
  for ( ; it != end; ++it)
    it->second->set_material(mat);
}


void
PhysicalModelInterface::add_submodel(const std::string& key, PhysicalModelInterface* pm)
{
  if (pm != NULL)
  {
    pm->set_simulator_id(get_simulator_id());
    pm->set_owner(get_owner());
    pm->set_material(get_material());
    pm->get_options().set_key(key);
    _submodels.insert(SubmodelMap::value_type(key, pm));
  }
}



void
PhysicalModelInterface::delete_submodel(const std::string& key)
{
  _submodels.erase(key);
}



//template <>
void
PhysicalModelInterface::_create_submodel(PhysicalModelInterface*& model,
    const std::string& type)
{
  model = NULL;

  string modname(type);

  // loop over all submodels
  ModelOptions::submodel_iterator it(get_options().submodels_begin(type));
  const ModelOptions::submodel_iterator end(get_options().submodels_end(type));

  for ( ; it != end; ++it)
  {
    if (model != NULL)
      throw ModelErrorException("Only one instance of submodel type \'"
          + type + "\' allowed");

    string modtype = ((it->second).get_option("type", (it->second).get_name()));
    (it->second).set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;

    // we try to create it from the same module
    model = create(modname, get_owner(), it->second, get_module_name());

    if (model == NULL)
    {
      // perhaps it uses 'type' or 'model' internally?
      model = create(it->first, get_owner(), it->second, get_module_name());
    }

    add_submodel(type, model);

    if (model == NULL)
    {
      ostringstream os;
      os << "Unknown physical model \'" << it->first << "\' (type \'"
        << modtype << "\')";
      throw InitFailedException(os.str());
    }
  }
}


//template <>
void
PhysicalModelInterface::_create_submodel(PhysicalModelInterface*& model,
    const std::string& type, const ModelOptions& default_opts)
{
  create_submodel(model, type);

  if (model == NULL)
  {
    string modname(type);

    ModelOptions opts(default_opts);
    string modtype = (opts.get_option("type", opts.get_name()));
    opts.set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;

    // we try to create it from the same module
    model = create(modname, get_owner(), opts, get_module_name());

    if (model == NULL)
    {
      // perhaps it uses 'type' or 'model' internally?
      model = create(type, get_owner(), opts, get_module_name());
    }

    add_submodel(type, model);

    if (model == NULL)
    {
      ostringstream os;
      os << "Unknown physical model \'" << type << "\' (type \'"
        << modtype << "\')";
      throw InitFailedException(os.str());
    }
  }
}


//template <>
void
PhysicalModelInterface::_create_submodels(std::vector<PhysicalModelInterface*>& models,
    const std::string& type)
{
  models.resize(0);

  // loop over all submodels
  ModelOptions::submodel_iterator it(get_options().submodels_begin(type));
  const ModelOptions::submodel_iterator end(get_options().submodels_end(type));

  for ( ; it != end; ++it)
  {
    string modname(type);

    string modtype = ((it->second).get_option("type", (it->second).get_name()));
    (it->second).set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;


    // we try to create it from the same module
    PhysicalModelInterface* mod =
        create(modname, get_owner(), it->second, get_module_name());

    if (mod == NULL)
    {
      // perhaps it uses 'type' or 'model' internally?
      mod = create(it->first, get_owner(), it->second, get_module_name());
    }

    add_submodel(type, mod);

    if (mod != NULL)
      models.push_back(mod);
    else
    {
      ostringstream os;
      os << "Unknown physical model \'" << it->first << "\' (type \'"
        << modtype << "\')";
      throw InitFailedException(os.str());
    }
  }
}



//template <>
void
PhysicalModelInterface::_create_submodels(std::vector<PhysicalModelInterface*>& models,
    const std::string& type, const ModelOptions& default_opts)
{
  create_submodels(models, type);

  if (models.empty())
  {
    PhysicalModelInterface* mod = NULL;

    string modname(type);

    ModelOptions opts(default_opts);
    string modtype = (opts.get_option("type", opts.get_name()));
    opts.set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;

    // we try to create it from the same module
    mod = create(modname, get_owner(), opts, get_module_name());

    if (mod == NULL)
    {
      // perhaps it uses 'type' or 'model' internally?
      mod = create(type, get_owner(), opts, get_module_name());
    }

    add_submodel(type, mod);

    if (mod != NULL)
      models.push_back(mod);
    else
    {
      ostringstream os;
      os << "Unknown physical model \'" << type << "\' (type \'"
        << modtype << "\')";
      throw InitFailedException(os.str());
    }
  }
}



void
PhysicalModelInterface::_create_submodels(void)
{
  // first call the user defined method
  prepare_submodels();
}


