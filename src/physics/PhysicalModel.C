// $Id$

#include "tibercad/base/tiber_config.h"
#include "tibercad/physics/PhysicalModel.h"
#include "tibercad/module/SimulationInterface.h"
#include "tibercad/physics/MaterialBoundary.h"
#include "tibercad/physics/Material.h"
#include "tibercad/io/Database.h"
#include "tibercad/math/TensorOperators.h"

#include "tibercad/utils/Utils.h"
#include "tibercad/physics/misc/Trap.h"
#include "tibercad/physics/misc/ParticleDensity.h"
#include "tibercad/physics/misc/PolarizationModel.h"

#include "tibercad/physics/semiconductormodels/ZbSemiconductor.h"
#include "tibercad/physics/semiconductormodels/WzSemiconductor.h"
#include "tibercad/physics/semiconductormodels/ZbDDsemiconductor.h"
#include "tibercad/physics/semiconductormodels/WzDDsemiconductor.h"

#include "tibercad/physics/semiconductormodels/SBZbCondBandBulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/SBWzCondBandBulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/SBuserHamiltonian.h"
#include "tibercad/physics/semiconductormodels/KPbulkHamiltonian.h"
#include "tibercad/physics/semiconductormodels/EFAbulkModel.h"




#include "tibercad/io/Messages.h"

#include <string>



using namespace std;


map<const string, ID>
PhysicalModel::_model_ids;

/*map< std::pair<const std::string,
               const std::string>,
     ID >
PhysicalModel::_unique_model_ids;*/


PhysicalModel::PhysicalModel(const ModelOptions& options)
  : TiberModelObject(options),
    _id(INVALID_ID),
    _simulator_id(INVALID_ID),
    _owner(nullptr),
    _bulk_material(nullptr),
    _database(nullptr),
    _module("")
{
}

PhysicalModel::~PhysicalModel(void)
{
  SubmodelIterator it(submodels_begin());
  const SubmodelIterator end(submodels_end());
  for ( ; it != end; ++it)
    destroy(it->second);
}



PhysicalModel*
PhysicalModel::_create(const string& name,
    const PhysicalObject* owner,
    const ModelOptions& options,
    const string& module)
{

  PhysicalModel* mod = nullptr;

  if (name == "semicond_zb")
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
  else if (name == "trap")
    mod = Trap::create(options);
  else if (name == "polarization")
    mod = PolarizationModel::create(options);
  else if (name == "particle_density")
    mod = ParticleDensity::create(options);


  if (mod == nullptr)
  {
    // first try without module directory
    if ((mod = create_from_library<PhysicalModel>(name, options, owner)) == 0)
    {
      if (module.size() >  0)
      {
        mod = create_from_library<PhysicalModel>(
            module + "/" + name, options, owner);
      }
    }
  }

  if (mod != nullptr)
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



PhysicalModel*
PhysicalModel::_create(create_t create_fnc, destroy_t destroy_fnc,
    const PhysicalObject* owner,
    const ModelOptions& options,
    const string& module)
{
  PhysicalModel* mod = dynamic_cast<PhysicalModel*>(
      create_from_function(create_fnc, destroy_fnc, options, owner));

  if (mod != nullptr)
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
PhysicalModel::_register_model(
    PhysicalModel* model)
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

  /*
  if (model->_has_unique_id) {
    const string mat_name = (model->get_material())->get_name();
    const string unique_name = model->get_name();

    unique_model_id_iterator it = _unique_model_ids.find(std::make_pair(unique_name, mat_name));

    if (it == _unique_model_ids.end())
    {
      id = _unique_model_ids.size() + 1;
      _unique_model_ids[std::make_pair(unique_name, mat_name)] = id;
      model->_unique_id = id;
    }
    else
      throw ModelErrorException("'" + name + "' models MUST have unique names. Name: '" + unique_name + "' is not unique in material '" + mat_name + "'" );
  }
  */

  model->_id = id;
}



void
PhysicalModel::set_owner(const PhysicalObject* owner)
{
  _owner = owner;
  //if ((_owner != nullptr) && (_owner->get_type() == PhysicalObject::BULK))
  //  _bulk_material = static_cast<const Material*>(_owner);
  //else
  //  _bulk_material = nullptr;
}



PhysicalModel*
PhysicalModel::copy(void) const
{
  PhysicalModel* new_copy = nullptr;

  // this is safe
  new_copy = this->create_new();

  if (new_copy != nullptr)
  {
    new_copy->_id = _id;
    new_copy->set_owner(_owner);
    new_copy->set_material(_bulk_material);
    new_copy->_database = _database;
    new_copy->set_name(get_name());
    new_copy->_set_type(get_type());
    new_copy->_simulator_id = _simulator_id;
    new_copy->_set_module_name(get_module_name());

    new_copy->copy_from(this);
  }

  return new_copy;
}



PhysicalModel*
PhysicalModel::create_new(void) const
{
  PhysicalModel* pmi = create_from_object(this, get_owner());
  if (pmi == nullptr)
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
PhysicalModel::get_default_name(void) const
{
  return Utils::extract_typename(typeid(*this));
}




const Database&
PhysicalModel::get_database(void)
{
  if (_database != nullptr)
    return(*_database);

  return _owner->get_database();
}

void
PhysicalModel::set_database(const Database* db)
{
  if (db != nullptr)
    _database = db;
  else
    _database = &_owner->get_database();
}


void
PhysicalModel::override_parameter_string(const std::string& name,
        std::string& s) const
{
  const PhysicalObject* mat = get_owner();

  if (mat != nullptr)
  {
    std::string code;

    SimulationInterface* sim =
      SimulationInterface::get_simulation(get_simulator_id());

    // first ask for the plain name
    s = mat->get_options().get_option(name, s);

    // first override
    if (sim != nullptr)
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


bool
PhysicalModel::override_material(void)
{
  bool overridden = false;

  if (get_options().has_submodel("override_material"))
  {
    ModelOptions opts =
      get_options().submodels_begin("override_material")->second;

    if (opts.get_option("copy_from_region", false))
    {
      opts = get_owner()->get_options();
      opts += get_options().submodels_begin("override_material")->second;
    }

    string mat = opts.get_name();
    Material* newmat = Material::create(mat, opts);

    get_options().delete_submodels("override_material");

    const PhysicalObject* obj = this->get_owner();
    newmat->add_model(this, this->get_simulator_id());
    //_database = new Database(get_database());
    _owner = obj;
    newmat->init();
    //_bulk_material = newmat;

    overridden = true;
  }

  return(overridden);
}




void
PhysicalModel::init(void)
{
  // first check if we want to override the material datafile

  if (!override_material())
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
    Messages::debug("init() of " + get_name() + " (" + get_material()->get_name() + ")");
    Messages m;
    m.indent();

    // Now we initialize all "official" submodels
    SubmodelIterator smit(submodels_begin());
    const SubmodelIterator smend(submodels_end());
    for ( ; smit != smend; ++smit)
    {
      PhysicalModel* pm = smit->second;

      //pm->_simulator_id = _simulator_id;
      // we pass overriden database, if present
      pm->_database = _database;
      Messages::debug("Initializing submodel " + smit->first + " ...");
      pm->init();
      Messages::debug("done (Initializing submodel " + smit->first + ")");
    }

    m.unindent();
    do_init();

    // dummy read
    get_option("crystal_structure", "");
    get_option("regions", "");
    get_option("type", "");
    get_option("name", "");

    get_options().check_unused(1);
    Messages::debug("finished init() of " + get_name() + " (" + get_material()->get_name() + ")");
  }
}




void
PhysicalModel::reinit(void)
{
  // reinit submodels
  SubmodelIterator smit(submodels_begin());
  const SubmodelIterator smend(submodels_end());
  for ( ; smit != smend; ++smit)
    smit->second->reinit();

  do_reinit();
}




void
PhysicalModel::reinit(const Elem* elem)
{
  // reinit submodels
  SubmodelIterator smit(submodels_begin());
  const SubmodelIterator smend(submodels_end());
  for ( ; smit != smend; ++smit)
    smit->second->reinit(elem);

  do_reinit(elem);
}




void
PhysicalModel::init_interface(const Material* comp_A,
    const Material* comp_B)
{
  // sometimes this is not a good idea, and it even may not make
  // sense at all in general
  //read_database();

  read_interface_database();

  // setup the submodels
  _create_submodels();
  Messages::debug("init_interface() of " + this->get_name() + " (" + get_owner()->get_name() + ")");

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
PhysicalModel::init_alloy(const PhysicalModel* comp_A,
    const PhysicalModel* comp_B, double xa)
{
  assert(typeid(*comp_A) == typeid(*comp_B));
  Messages::debug("init_alloy() of " + this->get_name() + " (" + get_owner()->get_name() + ")");

  //
  // NOTE:
  //   The current approach relies on the strong assumption that all models
  //   in all alloy components are ordered exactly the same way. This had
  //   better be changed in the future!

  // first check if we want to override the material datafile
  if (!override_material())
  {

    read_database();

    // some models might treat alloys in a special way
    // disable alloy mixing
    Database::AlloyMixing mixing = get_database().get_alloy_mixing();
    get_database().set_alloy_mixing(Database::NONE);
    read_database_alloy();
    get_database().set_alloy_mixing(mixing);


    // setup the submodels
    _create_submodels();

    Messages m;
    m.indent();

    SubmodelIterator it(submodels_begin());
    ConstSubmodelIterator it_A(comp_A->submodels_begin());
    ConstSubmodelIterator it_B(comp_B->submodels_begin());
    const SubmodelIterator end(submodels_end());
    for ( ; it != end; ++it, ++it_A, ++it_B)
    {
      PhysicalModel* pm = it->second;
      pm->init_alloy(it_A->second, it_B->second, xa);
    }
    m.unindent();

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
}


void
PhysicalModel::set_material(const Material* mat)
{
  _bulk_material = mat;

  SubmodelIterator it(_submodels.begin());
  const SubmodelIterator end(_submodels.end());
  for ( ; it != end; ++it)
    it->second->set_material(mat);
}


void
PhysicalModel::add_submodel(const std::string& key, PhysicalModel* pm)
{
  if (pm != nullptr)
  {
    pm->set_simulator_id(get_simulator_id());
    pm->set_owner(get_owner());
    pm->set_material(get_material());
    pm->_database = _database;
    pm->get_options().set_key(key);
    _submodels.insert(SubmodelMap::value_type(key, pm));
  }
}



void
PhysicalModel::delete_submodel(const std::string& key)
{
  _submodels.erase(key);
}


void
PhysicalModel::delete_submodel(SubmodelIterator it)
{
  _submodels.erase(it);
}

//template <>
void
PhysicalModel::_create_submodel(PhysicalModel*& model,
    const std::string& type)
{
  model = nullptr;

  string modname(type);

  // loop over all submodels
  ModelOptions::submodel_iterator it(get_options().submodels_begin(type));
  const ModelOptions::submodel_iterator end(get_options().submodels_end(type));

  for ( ; it != end; ++it)
  {
    if (model != nullptr)
      throw ModelErrorException("Only one instance of submodel type \'"
          + type + "\' allowed");

    string modtype = ((it->second).get_option("type", (it->second).get_name()));
    (it->second).set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;

    // we try to create it from the same module
    model = create(modname, get_owner(), it->second, get_module_name());

    if (model == nullptr)
    {
      // perhaps it uses 'type' or 'model' internally?
      model = create(it->first, get_owner(), it->second, get_module_name());
    }

    add_submodel(type, model);

    if (model == nullptr)
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
PhysicalModel::_create_submodel(PhysicalModel*& model,
    const std::string& type, const ModelOptions& default_opts)
{
  create_submodel(model, type);

  if (model == nullptr)
  {
    string modname(type);

    ModelOptions opts(default_opts);
    string modtype = (opts.get_option("type", opts.get_name()));
    opts.set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;

    // we try to create it from the same module
    model = create(modname, get_owner(), opts, get_module_name());

    if (model == nullptr)
    {
      // perhaps it uses 'type' or 'model' internally?
      model = create(type, get_owner(), opts, get_module_name());
    }

    add_submodel(type, model);

    if (model == nullptr)
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
PhysicalModel::_create_submodels(std::vector<PhysicalModel*>& models,
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
    PhysicalModel* mod =
        create(modname, get_owner(), it->second, get_module_name());

    if (mod == nullptr)
    {
      // perhaps it uses 'type' or 'model' internally?
      mod = create(it->first, get_owner(), it->second, get_module_name());
    }

    add_submodel(type, mod);

    if (mod != nullptr)
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
PhysicalModel::_create_submodels(std::vector<PhysicalModel*>& models,
    const std::string& type, const ModelOptions& default_opts)
{
  create_submodels(models, type);

  if (models.empty())
  {
    PhysicalModel* mod = nullptr;

    string modname(type);

    ModelOptions opts(default_opts);
    string modtype = (opts.get_option("type", opts.get_name()));
    opts.set_option("type", modtype);

    if (modtype.size() > 0)
      modname += string("_") + modtype;

    // we try to create it from the same module
    mod = create(modname, get_owner(), opts, get_module_name());

    if (mod == nullptr)
    {
      // perhaps it uses 'type' or 'model' internally?
      mod = create(type, get_owner(), opts, get_module_name());
    }

    add_submodel(type, mod);

    if (mod != nullptr)
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
PhysicalModel::_create_submodels(void)
{
  // first call the user defined method
  prepare_submodels();
}



void
PhysicalModel::alloy(Tensor4DSym& result, const Tensor4DSym& val_a,
    const Tensor4DSym& val_b, double xa)
{
 result = (1 - xa) * val_b + xa * val_a;
}

void
PhysicalModel::alloy(Tensor4DSym& result, const Tensor4DSym& val_a,
    const Tensor4DSym& val_b, double xa, const Tensor4DSym& bowing)
{
 result = (1 - xa) * val_b + xa * val_a - xa * (1 - xa) * bowing;
}


void
PhysicalModel::alloy(Tensor2& result, const Tensor2& val_a,
    const Tensor2& val_b, double xa)
{
  result = (1 - xa) * val_b + xa * val_a;
}

void
PhysicalModel::alloy(Tensor2& result, const Tensor2& val_a,
    const Tensor2& val_b, double xa, const Tensor2& bowing)
{
  result = (1 - xa) * val_b + xa * val_a - xa * (1 - xa) * bowing ;
}
