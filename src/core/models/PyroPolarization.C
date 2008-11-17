// $Id$

#include "PyroPolarization.h"
#include "Material.h"


PhysicalModelInterface*
PyroPolarization::create_new(void) const
{
  return new PyroPolarization();
}


PyroPolarization*
PyroPolarization::create(void)
{
  return new PyroPolarization();
}



PyroPolarization*
PyroPolarization::create(const Material* mat)
{
  if (mat == NULL)
    throw InitFailedException("Try to create PyroPolarization with invalid material.");

  PyroPolarization* pyro = NULL;

  if (mat->get_structure() == "wz")
    pyro = dynamic_cast<PyroPolarization*>(
        PhysicalModelInterface::create("pyropolarization_wz"));
  else
    pyro = dynamic_cast<PyroPolarization*>(
        PhysicalModelInterface::create("pyropolarization_zb"));
  
  return pyro;
}
