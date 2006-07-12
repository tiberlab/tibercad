// $Id$


#include "RecombinationModelInterface.h"
#include "SRHRecombination.h"
#include "DirectRecombination.h"

#include <string>


RecombinationModelInterface::RecombinationModelInterface(void)
  : DriftDiffusionModelInterface()
{
}

RecombinationModelInterface*
RecombinationModelInterface::create(const std::string& name)
{
  RecombinationModelInterface* rec = NULL;

  if (name == "SRH")
    rec = new SRHRecombination();
  else if (name == "direct")
    rec = new DirectRecombination();

  return rec;
}


RecombinationModelInterface*
RecombinationModelInterface::create(const std::string& name,
    const ModelOptions& options)
{
  RecombinationModelInterface* rec = create(name);

  if (rec != NULL)
    rec->set_model_options(options);

  return rec;
}


