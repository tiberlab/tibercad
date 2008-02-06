// $Id$

#include "tiber_config.h"
#include "PhysicalModelInterface.h"
#include "DLLoader.h"

#ifndef BUILD_TIBER_MODULES
#include "SRHRecombination.h"
#include "DirectRecombination.h"
#include "ExcitonGeneration.h"
#include "ExcitonDissociation.h"
#include "OpticalGeneration.h"

#include "ConstantMobility.h"
#include "DopingDependentMobility.h"
#include "FieldDependentMobility.h"

#include "SimpleSemiconductorModel.h"
#include "SemiconductorModel.h"
#include "StrainedSemiconductorModel.h"
#endif
#include "MobilityModelInterface.h"
#include "RecombinationModelInterface.h"

#include "ExcitonModel.h"

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


#include "ThermoelectricPower.h"
#include "ParticleThermalConductivity.h"

#include  "ZbLatticeThermalConductivity.h"
#include  "WzLatticeThermalConductivity.h"
#include  "HeatModel.h"
#include  "TightBindingModel.h"
#include  "PoissonModel.h"
#include  "ChargeDensityModel.h"
#include  "DielectricModel.h"

#include  "ZbOptDielectricConstant.h"
#include  "WzOptDielectricConstant.h"
#include  "MaxwellPhysicalModel.h"

#include "WzPyroPolarization.h"


#include <typeinfo>
#ifdef DEBUG
#include <iostream>
#endif
#include <string>


using namespace std;


map<const string, ID> 
PhysicalModelInterface::_model_ids;


PhysicalModelInterface*
PhysicalModelInterface::create(const string& name,
    const ModelOptions& options)
{

  PhysicalModelInterface* mod = NULL;

  // First we attempt to open a shared library
  //
  DLLoader::LibraryInterface iface;
  bool success = DLLoader::open_library(name, iface);

  create_t create_fnc = (create_t) iface.create_fnc;
  destroy_t destroy_fnc = (destroy_t) iface.destroy_fnc;

  if(success)
    // Try to create the object
    mod = create_fnc();
  else
  {
#ifndef BUILD_TIBER_MODULES
    if (name == "dd_rec_srh")
      mod = SRHRecombination::create();
    else if (name == "dd_rec_direct")
      mod = DirectRecombination::create();
    else if (name == "dd_rec_exciton_generation")
      mod = ExcitonGeneration::create();
    else if (name == "dd_rec_exciton_dissociation")
      mod = ExcitonDissociation::create();
    else if (name == "dd_rec_optical")
      mod = OpticalGeneration::create();
    else if (name == "dd_mob_constant")
      mod = ConstantMobility::create();
    else if (name == "dd_mob_doping_dependent")
      mod = DopingDependentMobility::create();
    else if (name == "dd_mob_field_dependent")
      mod = FieldDependentMobility::create();
    else if (name == "dd_simple")
      mod = SimpleSemiconductorModel::create();
    else if (name == "dd_unstrained")
      mod = SemiconductorModel::create();
    else if (name == "dd_strained")
      mod = StrainedSemiconductorModel::create();
    else if (name == "exmodel_simple")
#else    
    if (name == "exmodel_simple")
#endif
      mod = ExcitonModel::create();
    else if (name == "stiffness_zb")
      mod = ZbStiffness::create();
    else if (name == "stiffness_wz")
      mod = WzStiffness::create();
    else if (name == "piezo_zb")
      mod = ZbPiezoelectricity::create();
    else if (name == "piezo_wz")
      mod = WzPiezoelectricity::create();
    else if (name == "cryst_zb")
      mod = ZbRotatedCrystal::create();
    else if (name == "cryst_wz")
      mod = WzRotatedCrystal::create();
    else if (name == "macrostrain")
      mod = MacrostrainModel::create();
    else if (name == "semicond_zb")
      mod = ZbSemiconductor::create();
    else if (name == "semicond_wz")
      mod = WzSemiconductor::create();
    else if (name == "quantum_kp")
      mod = KPbulkHamiltonian::create();
    else if (name == "quantum_cond_band_zb")
      mod = SBZbCondBandBulkHamiltonian::create();
    else if (name == "quantum_cond_band_wz")
      mod = SBWzCondBandBulkHamiltonian::create();
    else if (name == "quantum_user")
      mod = SBuserHamiltonian::create();
    else if (name == "DDsemicond_zb")
      mod = ZbDDsemiconductor::create();
    else if (name == "DDsemicond_wz")
      mod = WzDDsemiconductor::create();
    else if (name == "EFAmodel")
      mod = EFAbulkModel::create();
    else if (name == "lat_therm_cond_zb")
      mod = ZbLatticeThermalConductivity::create();
    else if (name == "lat_therm_cond_wz")
      mod = WzLatticeThermalConductivity::create();
    else if (name == "thermal")
      mod = HeatModel::create();
    else if (name == "thermoelectric_power") 
      mod = ThermoelectricPower::create();
    else if  (name == "particle_thermal_conductivity")
      mod = ParticleThermalConductivity::create();
    else if  (name == "poisson")
      mod = PoissonModel::create();
    else if (name == "tightbinding")
       mod = TightBindingModel::create();
    else if  (name == "charge_density_model")
      mod = ChargeDensityModel::create();
    else if  (name == "dielectric_model")
      mod = DielectricModel::create();
    else if (name == "opt_dielectric_constant_zb")
      mod = ZbOptDielectricConstant::create();
    else if (name == "opt_dielectric_constant_wz")
      mod = WzOptDielectricConstant::create();
    else if (name == "maxwell")
      mod = MaxwellPhysicalModel::create();
    else if (name == "pyropolarization")
      mod = PyroPolarization::create();
    else if (name == "pyropolarization_wz")
      mod = WzPyroPolarization::create();
  }


  if (mod != NULL)
  {
    register_model(mod);

    mod->_libhandle = iface.handle;
    mod->_create = create_fnc;
    mod->_destroy = destroy_fnc;


    mod->set_options(options);

    //! set the name
    // 2007-08-17
    //    we don't set anymore a default name
    //string defaultname = mod->get_default_name();
    string defaultname = "";
    mod->_name = mod->_options.get_option("name", defaultname);
    mod->_options.delete_option("name");
#ifdef DEBUG
    cerr << "Add model (ID = " << mod->get_id() <<
      " name = " << mod->get_name() << " type_id = " <<
       mod->get_default_name() << ")\n";
#endif
  }

  return mod;
}


void
PhysicalModelInterface::destroy(PhysicalModelInterface* p)
{

  if (p != NULL)
  {
#ifdef DEBUG
    string id = Utils::extract_typename(typeid(*p));
    cerr << "Delete model (ID = " << p->get_id() <<
      " name = " << p->get_name() << " type_id = " << id << ")\n";
#endif

    libhandle_t libhandle = p->_libhandle;
    destroy_t destroy_fnc = p->_destroy;

    if (destroy_fnc != NULL)
      destroy_fnc(p);
    else
      delete p;

    if (libhandle != NULL)
      DLLoader::close_library(libhandle);
  }
}


template <typename T>
ID
PhysicalModelInterface::get_id_from_name(const string& name)
{
  ID id = 0;

  PhysicalModelInterface* rec = T::create(name);

  if (rec != NULL)
    id = rec->get_id();

  // rec is either valid or NULL
  destroy(rec);
  

  return id;
}


void
PhysicalModelInterface::register_model(
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
    new_copy->_material = _material;
    new_copy->_options = _options;
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



template <typename T>
T
PhysicalModelInterface::get_parameter(const std::string& name,
    T default_value) const
{
  T val(_options.get_option(name, default_value));

  const Material* mat = get_material();

  if (mat != NULL)
  {
    std::string code;

    SimulationInterface* sim = 
      SimulationInterface::get_simulation(get_simulator_id());

    // first ask for the plain name
    val = mat->get_options().get_option(name, val);

    // first override
    if (sim != NULL)
    {
      code = sim->get_name() + ".";
      val = mat->get_options().get_option(code + name, val);
    }

    // second override
    if (get_name() != "")
    {
      code += get_name() + ".";
      val = mat->get_options().get_option(code + name, val);
    }
  }

  return val;
}



template <typename T>
void
PhysicalModelInterface::get_parameter(const std::string& name,
    std::vector<T>& vec) const
{
  _options.get_option(name, vec);

  const Material* mat = get_material();

  if (mat != NULL)
  {
    std::string code;

    SimulationInterface* sim = 
      SimulationInterface::get_simulation(get_simulator_id());

    // first ask for the plain name
    mat->get_options().get_option(name, vec);

    // first override
    if (sim != NULL)
    {
      code = sim->get_name() + ".";
      mat->get_options().get_option(code + name, vec);
    }

    // second override
    if (get_name() != "")
    {
      code += get_name() + ".";
      mat->get_options().get_option(code + name, vec);
    }
  }
}




// explicit instantiations

template ID
PhysicalModelInterface::get_id_from_name<RecombinationModelInterface>(
    const string& name);

template ID
PhysicalModelInterface::get_id_from_name<MobilityModelInterface>(
    const string& name);







template
double
PhysicalModelInterface::get_parameter<double>(const string& name,
    double val) const;

template
int
PhysicalModelInterface::get_parameter<int>(const string& name,
    int val) const;

template
unsigned int
PhysicalModelInterface::get_parameter<unsigned int>(const string& name,
    unsigned int val) const;

template
short
PhysicalModelInterface::get_parameter<short>(const string& name,
    short val) const;


template
bool
PhysicalModelInterface::get_parameter<bool>(const string& name,
    bool val) const;

template
char
PhysicalModelInterface::get_parameter<char>(const string& name,
    char val) const;

template
string
PhysicalModelInterface::get_parameter<string>(const string& name,
    string val) const;

template
const char*
PhysicalModelInterface::get_parameter<const char*>(const string& name,
    const char* val) const;







template
void
PhysicalModelInterface::get_parameter<double>(const string& name,
    vector<double>& vec) const;

template
void
PhysicalModelInterface::get_parameter<int>(const string& name,
    vector<int>& vec) const;

template
void
PhysicalModelInterface::get_parameter<unsigned int>(const string& name,
    vector<unsigned int>& vec) const;

template
void
PhysicalModelInterface::get_parameter<short>(const string& name,
    vector<short>& vec) const;


template
void
PhysicalModelInterface::get_parameter<bool>(const string& name,
    vector<bool>& vec) const;

template
void
PhysicalModelInterface::get_parameter<char>(const string& name,
    vector<char>& vec) const;

template
void
PhysicalModelInterface::get_parameter<string>(const string& name,
    vector<string>& vec) const;


