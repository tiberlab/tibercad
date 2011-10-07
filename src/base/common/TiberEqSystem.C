// $Id$

#include "TiberEqSystem.h"
#include "TiberNonlinearSystem.h"
#include "TiberLinearSystem.h"
#include "InitFailedException.h"

#include <numeric_vector.h>


TiberEqSystem::TiberEqSystem(void)
  : _type(LINEAR),
    _l2_weight(NULL),
    _linfty_weight(NULL)
{
}



TiberEqSystem*
TiberEqSystem::create(EquationSystems& es,
    const std::string& sysname, const std::string& type,
    const ModelOptions& options)
{
  SystemType sys_type = LINEAR;

  if (type == "linear")
    sys_type = LINEAR;
  else if (type == "nonlinear")
    sys_type = NONLINEAR;
  else
    throw InitFailedException("Unknown equation system type " + sys_type);

  return create(es, sysname, sys_type, options);
}



TiberEqSystem*
TiberEqSystem::create(EquationSystems& es,
    const std::string& sysname, SystemType type,
    const ModelOptions& options)
{
  TiberEqSystem* sys = NULL;

  switch (type)
  {
    case LINEAR:
      sys = TiberLinearSystem::create(es, sysname, options);
      break;

    case NONLINEAR:
      sys = TiberNonlinearSystem::create(es, sysname, options);
      break;

    default:
      throw InitFailedException("Unknown equation system type");
  }

  return sys;
}
  


void
TiberEqSystem::set_weight(const NumericVector<double>* weight, NormType norm)
{
  switch (norm)
  {
    case l2_NORM:
      _l2_weight = weight;
      break;

    case MAX_NORM:
      _linfty_weight = weight;
      break;

    default:
      break;
  }
}


double
TiberEqSystem::calculate_norm(NumericVector<double>* vec, NormType norm)
{
  double result = 0;

  switch (norm)
  {
    case l2_NORM:
      if (_l2_weight != NULL)
        vec->pointwise_mult(*vec, *_l2_weight);
      result = vec->l2_norm();
      break;

    case MAX_NORM:
      if (_linfty_weight != NULL)
        vec->pointwise_mult(*vec, *_linfty_weight);
      result = vec->linfty_norm();
      break;

    default:
      break;
  }

  return result;
}


const NumericVector<double>*
TiberEqSystem::get_weight(NormType norm) const
{
  switch (norm)
  {
    case l2_NORM:
      return _l2_weight;
      break;

    case MAX_NORM:
      return _linfty_weight;
      break;

    default:
      break;
  }
  return NULL;
}
