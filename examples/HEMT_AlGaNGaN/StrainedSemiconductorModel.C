#include "StrainedSemiconductorModel.h"

#include "elem.h"
#include "macrostrain.h"
#include "tensor.h"

#include <iostream>

using namespace DriftDiffusionDefs;

StrainedSemiconductorModel::StrainedSemiconductorModel(
    Macrostrain* strain)
  : SimpleSemiconductorModel(),
    _strain(strain)
{
}


StrainedSemiconductorModel::StrainedSemiconductorModel(
    const StrainedSemiconductorModel& model)
  : SimpleSemiconductorModel(model),
    _strain(model._strain)
{
}


void
StrainedSemiconductorModel::calculate_all(
    double potential, double fermi_e, double fermi_h,
    const Point& p, const Elem* elem,
    int coupling)
{
  static const Elem* elem_old = NULL;

  if ((elem_old != elem) && (elem != NULL))
  {
    assert(elem != NULL);
    Tensor1 pol = _strain->get_piezopolarization(elem);
    polarization(0) = pol(1);
    polarization(1) = pol(2);
    polarization(2) = pol(3);
    elem_old = elem;
  }

  switch (coupling & BOTH)
  {
    case ELECTRONS:
      SimpleSemiconductorModel::calculate_all<ELECTRONS>(
          potential, fermi_e, fermi_h);
      break;
    case HOLES:
      SimpleSemiconductorModel::calculate_all<HOLES>(
          potential, fermi_e, fermi_h);
      break;
    default:
      SimpleSemiconductorModel::calculate_all<BOTH>(
          potential, fermi_e, fermi_h);
      break;
  }
}
