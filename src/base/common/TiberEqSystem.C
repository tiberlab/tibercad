// $Id$

#include "TiberEqSystem.h"
#include "TiberNonlinearSystem.h"
#include "TiberLinearSystem.h"
#include "InitFailedException.h"

#include <numeric_vector.h>
#include <dense_vector.h>
#include <dense_matrix.h>

#include <cassert>


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



System*
TiberEqSystem::get_libmesh_system(void)
{
  System* sys = NULL;
  switch (get_type())
  {
    case LINEAR:
      sys = static_cast<TiberLinearSystem*>(this);
      break;

    case NONLINEAR:
      sys = static_cast<TiberNonlinearSystem*>(this);
      break;
  }

  return sys;
}



void
TiberEqSystem::exclude_dofs(DenseMatrix<double>& mat,
    const std::vector<unsigned int>& dof_indices, const Elem* elem)
{

  // if there are no excluded Dofs we have nothing to do
  if (!_excluded_dofs.empty())
  {
    assert(dof_indices.size() == mat.m());
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
    {
      if (_excluded_dofs.count(dof_indices[i]))
      {
        for (unsigned int j = 0; j < mat.n(); ++j)
          mat(i,j) = 0;

        mat(i,i) = 1;
      }
      else if ((elem != NULL) && _excluded_region_ids.count(elem->subdomain_id()) &&
          _interface_dofs.count(dof_indices[i]))
      {
        // this DOF is on an interface, so it must not be messed up
        mat(i,i) = 0;
      }
    }
  }
}



void
TiberEqSystem::exclude_dofs(DenseVector<double>& vec,
    const std::vector<unsigned int>& dof_indices, const Elem* elem)
{
  // if there are no excluded Dofs we have nothing to do
  if (!_excluded_dofs.empty())
  {
    assert(dof_indices.size() == vec.size());
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
    {
      if (_excluded_dofs.count(dof_indices[i]))
        vec(i) = 0;
      else if ((elem != NULL) && _excluded_region_ids.count(elem->subdomain_id()) &&
          _interface_dofs.count(dof_indices[i]))
      {
        // this DOF is on an interface, so it must not be messed up
        vec(i) = 0;
      }
    }
  }
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
  double old;
  switch (norm)
  {
    case l2_NORM:
      old = vec->l2_norm();
      if (_l2_weight != NULL)
        vec->pointwise_mult(*vec, *_l2_weight);
      result = vec->l2_norm();
      //std::cerr << "l2 ";
      break;

    case MAX_NORM:
      old = vec->linfty_norm();
      if (_linfty_weight != NULL)
        vec->pointwise_mult(*vec, *_linfty_weight);
      result = vec->linfty_norm();
      //std::cerr << "linfty ";
      break;

    default:
      break;
  }
      //std::cerr << old << "  " << result << std::endl;

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
