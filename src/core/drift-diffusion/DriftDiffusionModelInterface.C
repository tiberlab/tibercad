// $Id$

#include "DriftDiffusionModelInterface.h"

#include <typeinfo>
#include <iostream>
#include <string>

std::map<const std::string, ID> 
DriftDiffusionModelInterface::_model_ids;

void
DriftDiffusionModelInterface::register_model(void)
{
  const std::string name = typeid(*this).name();
  model_id_iterator it = _model_ids.find(name);

  // set the model ID (create a new one if the model didn't exist yet)
  if (it == _model_ids.end())
  {
    _id = _model_ids.size() + 1;
    _model_ids[name] = _id;
  }
  else
    _id = it->second;

  // FIXME remove this
  std::cerr << name << " created. ID = " << _id << "\n";
}

