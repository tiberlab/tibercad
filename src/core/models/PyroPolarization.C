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



void
PyroPolarization::copy_from(const PhysicalModelInterface *rhs)
{
  const PyroPolarization* temp = dynamic_cast<const PyroPolarization*>(rhs);

  _polarization = temp->_polarization;
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
    pyro = new PyroPolarization();

  return pyro;
}
