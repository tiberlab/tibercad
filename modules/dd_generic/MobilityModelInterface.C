// $Id$

#include "MobilityModelInterface.h"
#include "tibercad/physics/Material.h"


MobilityModelInterface*
MobilityModelInterface::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  return PhysicalModel::create<MobilityModelInterface>("mobility_" + name, mat, options);
}
