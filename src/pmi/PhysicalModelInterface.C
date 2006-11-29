// $Id$

#include "PhysicalModelInterface.h"

// the models
// this will be done in a more elegant way
#include "SRHRecombination.h"
#include "DirectRecombination.h"
#include "ExcitonGeneration.h"
#include "ExcitonDissociation.h"
#include "OpticalGeneration.h"

#include "ConstantMobility.h"

#include "ExcitonModel.h"

#include "Utils.h"

#include "SimpleSemiconductorModel.h"
#include "SemiconductorModel.h"
#include "StrainedSemiconductorModel.h"

#include <typeinfo>
#include <iostream>
#include <string>

std::map<const std::string, ID> 
PhysicalModelInterface::_model_ids;


PhysicalModelInterface*
PhysicalModelInterface::create(const std::string& name,
    const ModelOptions& options)
{
  PhysicalModelInterface* mod = NULL;

  if (name == "rec_SRH")
    mod = SRHRecombination::create();
  else if (name == "rec_direct")
    mod = DirectRecombination::create();
  else if (name == "rec_exciton_generation")
    mod = ExcitonGeneration::create();
  else if (name == "rec_exciton_dissociation")
    mod = ExcitonDissociation::create();
  else if (name == "rec_optical")
    mod = OpticalGeneration::create();
  else if (name == "mob_constant")
    mod = ConstantMobility::create();
  else if (name == "ddmodel_simple")
    mod = SimpleSemiconductorModel::create();
  else if (name == "ddmodel_unstrained")
    mod = SemiconductorModel::create();
  else if (name == "ddmodel_strained")
    mod = StrainedSemiconductorModel::create();
  else if (name == "exmodel_simple")
    mod = ExcitonModel::create();

  register_model(mod);

  if (mod != NULL)
  {
    mod->set_options(options);

    //! set the name
    std::string defaultname = Utils::extract_typename(typeid(*mod));
    mod->_name = mod->_options.get_option("name", defaultname);
    mod->_options.delete_option("name");

    std::cerr << "Added model (ID = " << mod->get_id() <<
      " name = " << mod->get_name() << ")\n";
  }

  return mod;
}

template <typename T>
ID
PhysicalModelInterface::get_id_from_name(const std::string& name)
{
  ID id = 0;

  PhysicalModelInterface* rec = T::create(name);

  if (rec != NULL)
    id = rec->get_id();

  // rec is either valid or NULL
  delete rec;

  return id;
}


void
PhysicalModelInterface::register_model(
    PhysicalModelInterface* model)
{
  if (model != NULL)
  {
    const std::string name = typeid(*model).name();
    model_id_iterator it = _model_ids.find(name);

    ID id;

    // set the model ID (create a new one if the model didn't exist yet)
    if (it == _model_ids.end())
    {
      id = _model_ids.size() + 1;
      _model_ids[name] = id;

      // FIXME remove this
      std::cerr << "registered new ";
    }
    else
      id = it->second;

    model->_id = id;

    // FIXME remove this
    std::cerr << name << " model. ID = " << id << "\n";
  }
}

PhysicalModelInterface*
PhysicalModelInterface::copy(void) const
{
  PhysicalModelInterface* new_copy = NULL;

  new_copy = this->create_new();
  const std::string name = typeid(*new_copy).name();
  std::cerr << name << "\n";

  if (new_copy != NULL)
    new_copy->copy_from(this);

  return new_copy;
}


// explict instantiations

template ID
PhysicalModelInterface::get_id_from_name<RecombinationModelInterface>(
    const std::string& name);

template ID
PhysicalModelInterface::get_id_from_name<MobilityModelInterface>(
    const std::string& name);
