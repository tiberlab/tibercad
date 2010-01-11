// $Id$

#include "PyroPolarization.h"
#include "Material.h"


PhysicalModelInterface*
PyroPolarization::create_new(void) const
{
  return new PyroPolarization(get_options());
}


PyroPolarization*
PyroPolarization::create(const ModelOptions& options)
{
  return new PyroPolarization(options);
}



PyroPolarization*
PyroPolarization::create(const Material* mat, const ModelOptions& options)
{
  if (mat == NULL)
    throw InitFailedException("Try to create PyroPolarization with invalid material.");

  PyroPolarization* pyro = NULL;

  if (mat->get_structure() == "wz")
    pyro = dynamic_cast<PyroPolarization*>(
        PhysicalModelInterface::create("pyropolarization_wz", options));
  else
    pyro = dynamic_cast<PyroPolarization*>(
        PhysicalModelInterface::create("pyropolarization_zb", options));
  
  return pyro;
}
