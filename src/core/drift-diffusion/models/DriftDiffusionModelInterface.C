// $Id$

#include "DriftDiffusionModelInterface.h"

// the models
// this will be done in a more elegant way
#include "SRHRecombination.h"
#include "DirectRecombination.h"
#include "ExcitonGeneration.h"
#include "ExcitonDissociation.h"
//#include "OpticalGeneration.h"

#include <typeinfo>
#include <iostream>
#include <string>

std::map<const std::string, ID> 
DriftDiffusionModelInterface::_model_ids;

DriftDiffusionModelInterface*
DriftDiffusionModelInterface::create(const std::string& name)
{
  DriftDiffusionModelInterface* rec = NULL;

  if (name == "SRH")
    rec = new SRHRecombination();
  else if (name == "direct")
    rec = new DirectRecombination();
  else if (name == "exciton_generation")
    rec = new ExcitonGeneration();
  else if (name == "exciton_dissociation")
    rec = new ExcitonDissociation();
  //else if (name == "optical")
  //  rec = new OpticalGeneration();

  register_model(rec);

  return rec;
}


DriftDiffusionModelInterface*
DriftDiffusionModelInterface::create(const std::string& name,
    const ModelOptions& options)
{
  DriftDiffusionModelInterface* rec = create(name);

  if (rec != NULL)
    rec->set_model_options(options);

  return rec;
}

ID
DriftDiffusionModelInterface::get_id(const std::string& name)
{
  ID id = 0;

  DriftDiffusionModelInterface* rec = create(name);

  if (rec != NULL)
    id = rec->get_id();

  // rec is either valid or NULL
  delete rec;

  return id;
}


void
DriftDiffusionModelInterface::register_model(
    DriftDiffusionModelInterface* model)
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
    }
    else
      id = it->second;

    model->_id = id;

    // FIXME remove this
    std::cerr << name << " registered. ID = " << id << "\n";
  }
}


