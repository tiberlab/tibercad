#include "StrainedSemiconductorModel.h"

#include "elem.h"
#include "macrostrain.h"
#include "tensor.h"

#include <iostream>

using namespace DriftDiffusionDefs;

StrainedSemiconductorModel::StrainedSemiconductorModel(
    Macrostrain* strain)
  : SimpleSemiconductorModel(),
    _strain(strain),
    _ignore_strain(false)
{
}


StrainedSemiconductorModel::StrainedSemiconductorModel(
    const StrainedSemiconductorModel& model)
  : SimpleSemiconductorModel(model),
    _strain(model._strain),
    _ignore_strain(model._ignore_strain)
{
}

void
StrainedSemiconductorModel::prepare_element_data(void)
{
  static const Elem* elem_old = NULL;

  if (!_ignore_strain)
  {
    if ((elem_old != elem) && (elem != NULL))
    {
      Tensor1 pol = _strain->get_piezopolarization(elem);
      polarization(0) = pol(1);
      polarization(1) = pol(2);
      polarization(2) = pol(3);
      elem_old = elem;
    }
  }

  SimpleSemiconductorModel::prepare_element_data();
}

