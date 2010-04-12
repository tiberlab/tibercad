// $Id$

#include "tiber_config.h"
#include "PhysicalModelInterface.h"
#include "Material.h"
#include "Variable.h"

#ifndef BUILD_TIBER_MODULES
#include "SRHRecombination.h"
#include "AugerRecombination.h"
#include "DirectRecombination.h"
#include "ExcitonGeneration.h"
#include "ExcitonDissociation.h"
#include "OpticalGeneration.h"

#include "ConstantMobility.h"
#include "DopingDependentMobility.h"
#include "FieldDependentMobility.h"
#include "FieldAssistedMobility.h"

#include "SimpleSemiconductorModel.h"
#include "SemiconductorModel.h"

#include "ExcitonModel.h"

#include "DSSCModel.h"
#endif


#include "Utils.h"

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



#include  "ZbLatticeThermalConductivity.h"
#include  "WzLatticeThermalConductivity.h"
#include  "HeatModel.h"
#include  "DftbModel.h"
#include  "EtbModel.h"
#include  "PoissonModel.h"
#include  "ChargeDensityModel.h"
#include  "DielectricModel.h"
#include  "HeatSourceInterface.h"
#include  "DriftDiffusionHeatSource.h"
#include  "ConstantHeatSource.h"

#include  "ZbOptDielectricConstant.h"
#include  "WzOptDielectricConstant.h"
#include  "MaxwellPhysicalModel.h"

#include "WzPyroPolarization.h"

#include "PhononModel.h"
#include "ZbFreeDynamicalMatrix.h"
#include "ZbStrainDynamicalMatrix.h"
#include "ZbRamanTensor.h"

#include <ZbPiezoelectricModel.h>
#include <WzPiezoelectricModel.h>


#include "Messages.h"

#include <typeinfo>
#include <string>


using namespace std;


map<const string, ID>
PhysicalModelInterface::_model_ids;



PhysicalModelInterface::~PhysicalModelInterface(void)
{
  string id = Utils::extract_typename(typeid(*this));
  ostringstream os;
  os << "Delete model (ID = " << get_id() <<
    " name = " << get_name() << " type = " << get_type() << ")";
  Messages::debug(os.str());

  SubmodelIterator it(submodels_begin());
  const SubmodelIterator end(submodels_end());
  for ( ; it != end; ++it)
    destroy(it->second);
}



PhysicalModelInterface*
PhysicalModelInterface::create(const string& name,
    const ModelOptions& options, const string& module)
{

  // NOTE: for bulk models options contains the crystal structure

  PhysicalModelInterface* mod = NULL;

#ifndef BUILD_TIBER_MODULES
  if (name == "recombination_srh")
    mod = SRHRecombination::create(options);
  else if (name == "recombination_auger")
    mod = AugerRecombination::create(options);
  else if (name == "recombination_direct")
    mod = DirectRecombination::create(options);
  else if (name == "recombination_exciton_generation")
    mod = ExcitonGeneration::create(options);
  else if (name == "recombination_exciton_dissociation")
    mod = ExcitonDissociation::create(options);
  else if (name == "generation_optical")
    mod = OpticalGeneration::create(options);
  else if (name == "mobility_constant")
    mod = ConstantMobility::create(options);
  else if (name == "mobility_doping_dependent")
    mod = DopingDependentMobility::create(options);
  else if (name == "mobility_field_dependent")
    mod = FieldDependentMobility::create(options);
  else if (name == "mobility_field_assisted")
    mod = FieldAssistedMobility::create(options);
  else if (name == "driftdiffusion_simple")
    mod = SimpleSemiconductorModel::create(options);
  else if (name == "driftdiffusion_default")
    mod = SemiconductorModel::create(options);
  else if (name == "dssc_default")
    mod = DSSCModel::create(options);
  else if (name == "ex_simple")
    mod = ExcitonModel::create(options);
  else if (name == "stiffness_zb")
#else
  if (name == "stiffness_zb")
#endif
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
  else if (name == "lat_therm_cond_zb")
    mod = ZbLatticeThermalConductivity::create(options);
  else if (name == "lat_therm_cond_wz")
    mod = WzLatticeThermalConductivity::create(options);
  else if (name == "thermal")
    mod = HeatModel::create(options);
  else if  (name == "poisson")
    mod = PoissonModel::create(options);
  else if (name == "dftb")
    mod = DftbModel::create(options);
  else if (name == "etb")
    mod = ETBModel::create(options);
  else if  (name == "charge_density_model")
    mod = ChargeDensityModel::create(options);
  else if  (name == "dielectric_model")
    mod = DielectricModel::create(options);
  else if (name == "opt_dielectric_constant_zb")
    mod = ZbOptDielectricConstant::create(options);
  else if (name == "opt_dielectric_constant_wz")
    mod = WzOptDielectricConstant::create(options);
  else if (name == "maxwell")
    mod = MaxwellPhysicalModel::create(options);
  else if (name == "pyropolarization_zb")
    mod = PyroPolarization::create(options);
  else if (name == "pyropolarization_wz")
    mod = WzPyroPolarization::create(options);
  else if  (name == "drift_diffusion_dissipation")
    mod = DriftDiffusionHeatSource::create(options);
  else if  (name == "phonon")
    mod = PhononModel::create(options);
  else if  (name == "constant_heat_source")
    mod = ConstantHeatSource::create(options);
  else if  (name == "free_dynamical_matrix_zb")
    mod = ZbFreeDynamicalMatrix::create(options);
  else if  (name == "strain_dependent_zb")
    mod = ZbStrainDynamicalMatrix::create(options);
  else if  (name == "raman_tensor_zb")
    mod = ZbRamanTensor::create(options);
  else if  (name == "piezoelectric_model_zb")
    mod = ZbPiezoelectricModel::create(options);
  else if  (name == "piezoelectric_model_wz")
    mod = WzPiezoelectricModel::create(options);


  if (mod == NULL)
  {
    // first try in the module directory
    if ((module.size() > 0) && ((mod = create_from_library<PhysicalModelInterface>(
        module + "/" + name, options)) == 0))
    {
      mod = create_from_library<PhysicalModelInterface>(name, options);
    }
  }

  if (mod != NULL)
  {
    _register_model(mod);

    // we let it know what's its identifier
    mod->_set_type(name);

    mod->_set_module_name(module);

    //! set the name
    // 2007-08-17
    //    we don't set anymore a default name
    //string defaultname = mod->get_default_name();
    string defaultname = "";
    mod->_name = mod->get_options().get_option("name", defaultname);
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
PhysicalModelInterface::create(create_t create_fnc, destroy_t destroy_fnc,
    const ModelOptions& options)
{
  PhysicalModelInterface* mod = dynamic_cast<PhysicalModelInterface*>(
      create_from_function(create_fnc, destroy_fnc, options));

  if (mod != NULL)
  {
    _register_model(mod);

    //! set the name
    // 2007-08-17
    //    we don't set anymore a default name
    //string defaultname = mod->get_default_name();
    string defaultname = "";
    mod->_name = mod->get_options().get_option("name", defaultname);
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




PhysicalModelInterface*
PhysicalModelInterface::copy(void) const
{
  PhysicalModelInterface* new_copy = NULL;

  new_copy = this->create_new();

  if (new_copy != NULL)
  {
    new_copy->_id = _id;
    new_copy->_owner = _owner;
    new_copy->_name = _name;
    new_copy->_simulator_id = _simulator_id;

    new_copy->copy_from(this);
  }

  return new_copy;
}


string
PhysicalModelInterface::get_default_name(void) const
{
  return Utils::extract_typename(typeid(*this));
}




Database&
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




Material*
PhysicalModelInterface::get_material(void)
{
  Material* mat = NULL;
  if (_owner->get_type() == PhysicalObject::BULK)
    mat = static_cast<Material*>(_owner);
  return mat;
}



const Material*
PhysicalModelInterface::get_material(void) const
{
  const Material* mat = NULL;
  if (_owner->get_type() == PhysicalObject::BULK)
    mat = static_cast<const Material*>(_owner);
  return mat;
}



void
PhysicalModelInterface::init(void)
{
  read_database();

  _create_submodels();

  // Now we initialize all "official" submodels
  SubmodelIterator smit(submodels_begin());
  const SubmodelIterator smend(submodels_end());
  for ( ; smit != smend; ++smit)
  {
    PhysicalModelInterface* pm = smit->second;
    pm->_simulator_id = _simulator_id;
    Messages::debug("Initializing " + smit->first + " ...");
    pm->init();
  }

  do_init();
}


void
PhysicalModelInterface::init_interface(const Material* comp_A,
    const Material* comp_B)
{

  read_database();

  read_interface_database();

  // setup the submodels
  _create_submodels();

  SubmodelIterator it(submodels_begin());
  const SubmodelIterator end(submodels_end());

  for ( ; it != end; ++it)
    it->second->init_interface(comp_A, comp_B);

  // some models might treat alloys in a special way
  do_init_interface(comp_A, comp_B);

  do_init();
}


void
PhysicalModelInterface::init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  assert(typeid(*comp_A) == typeid(*comp_B));


  read_database();

  // some models might treat alloys in a special way
  // disable alloy mixing
  Database::AlloyMixing mixing = get_database().get_alloy_mixing();
  get_database().set_alloy_mixing(Database::NONE);
  read_database_alloy();
  get_database().set_alloy_mixing(mixing);


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
}



void
PhysicalModelInterface::add_submodel(const std::string& key, PhysicalModelInterface* pm)
{
  assert(pm != NULL);
  pm->set_simulator_id(get_simulator_id());
  pm->set_owner(get_owner());
  _submodels.insert(SubmodelMap::value_type(key, pm));
}



void
PhysicalModelInterface::_create_submodels(void)
{
  // first call the user defined method
  create_submodels();

  // loop over all submodels
  ModelOptions::submodel_iterator it(get_options().submodels_begin());
  const ModelOptions::submodel_iterator end(get_options().submodels_end());

  while (it != end)
  {
    string name(it->first);
    string type((it->second).get_option("model", ""));
    type = ((it->second).get_option("type", type));
    if (type.size() > 0)
      name += string("_") + type;

    // we try to create it from the same module
    PhysicalModelInterface* pm = create(name, it->second, get_module_name());

    if (pm == NULL)
    {
      ostringstream os;
      os << "Unknown physical model \'" << it->first << "\' (type \'"
        << type << "\')";
      throw InitFailedException(os.str());
    }

    add_submodel(it->first, pm);

    // a temporary iterator as we cannot delete the loop iterator
    ModelOptions::submodel_iterator tmp_it(it);

    // next entry
    ++it;

    // we delete the options from the ModelOptions object
    get_options().delete_submodel(tmp_it);
  }
}


