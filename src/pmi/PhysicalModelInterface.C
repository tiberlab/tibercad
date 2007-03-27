// $Id$

#include "PhysicalModelInterface.h"
#include "DLLoader.h"

// the models
// this will be done in a more elegant way
#include "SRHRecombination.h"
#include "DirectRecombination.h"
#include "ExcitonGeneration.h"
#include "ExcitonDissociation.h"
#include "OpticalGeneration.h"

#include "ConstantMobility.h"
#include "DopingDependentMobility.h"

#include "ExcitonModel.h"

#include "Utils.h"

#include "SimpleSemiconductorModel.h"
#include "SemiconductorModel.h"
#include "StrainedSemiconductorModel.h"

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


#include  "ZbLatticeThermalConductivity.h"
#include  "WzLatticeThermalConductivity.h"
#include  "HeatModel.h"

#include <dlfcn.h>
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

  if(success)
  {
    create_t create_fnc = (create_t) iface.create_fnc;
    destroy_t destroy_fnc = (destroy_t) iface.destroy_fnc;

    // Try to create the object
    mod = create_fnc();

    mod->_libhandle = iface.handle;
    mod->_create = create_fnc;
    mod->_destroy = destroy_fnc;
  }
  else
  {
    if (name == "dd_rec_SRH")
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
    else if (name == "ddmodel_simple")
      mod = SimpleSemiconductorModel::create();
    else if (name == "ddmodel_unstrained")
      mod = SemiconductorModel::create();
    else if (name == "ddmodel_strained")
      mod = StrainedSemiconductorModel::create();
    else if (name == "exmodel_simple")
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
  }


  if (mod != NULL)
  {
    register_model(mod);

    mod->set_options(options);

    //! set the name
    string defaultname = mod->get_default_name();
    mod->_name = mod->_options.get_option("name", defaultname);
    mod->_options.delete_option("name");
#ifdef DEBUG
    cout << "Add model (ID = " << mod->get_id() <<
      " name = " << mod->get_name() << " type_id = " <<
      defaultname << ")\n";
#endif
  }

  return mod;
}


void
PhysicalModelInterface::destroy(PhysicalModelInterface* p)
{
  // TODO call destroy of the model module

  if (p != NULL)
  {
#ifdef DEBUG
    string id = Utils::extract_typename(typeid(*p));
    cout << "Delete model (ID = " << p->get_id() <<
      " name = " << p->get_name() << " type_id = " << id << ")\n";
#endif
    libhandle_t libhandle = p->_libhandle;
    destroy_t destroy_fnc = p->_destroy;

    if (destroy_fnc != NULL)
      destroy_fnc(p);
    else
      delete p;

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




// explicit instantiations

template ID
PhysicalModelInterface::get_id_from_name<RecombinationModelInterface>(
    const string& name);

template ID
PhysicalModelInterface::get_id_from_name<MobilityModelInterface>(
    const string& name);
