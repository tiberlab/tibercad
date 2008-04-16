#include "StrainSimulation.h"
#include "MacrostrainModelInterface.h"
#include "ModelOptions.h"

PhysicalModel*
StrainSimulation::create_physical_model(const ModelOptions& options,
    const Material* mat) const throw (ModelErrorException)
{
  const std::string& modelname = options.get_option("model", "macrostrain");

  MacrostrainModelInterface* model =
    MacrostrainModelInterface::create(modelname, options);

  if (model == NULL)
    throw ModelErrorException(
        "StrainSimulation: No such physical model: " + modelname);

  return model;
}
