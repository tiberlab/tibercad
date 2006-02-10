// $Id: BoundaryDescriptor.C 14 2006-01-27 11:44:23Z maufder $

#include "BoundaryDescriptor.h"

const double
BoundaryDescriptor::_penalty_value = 1e6;


const std::vector<double>
BoundaryDescriptor::get_scaled_to_normal_derivative(
    const std::string& variable) const
{
  std::vector<double> coeff(2);

  const_iterator it = find(variable);
  const std::vector<double>& orig = it->second;

  if (std::fabs(orig[1]) < 1 / _penalty_value)
  {
    coeff[0] = orig[0] * _penalty_value;
    coeff[1] = orig[2] * _penalty_value;
  }

  return coeff;
}


