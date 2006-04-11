// $Id$

#include "BoundaryDescriptor.h"

#include <iostream>

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

void
BoundaryDescriptor::print_info() const
{
  const_iterator it = begin();
  const const_iterator it_end = end();

  for ( ; it != it_end; ++it)
  {

    const std::vector<double>& coeff = it->second;
    std::cout << it->first << ": "
      << " a = " << coeff[0]
      << " b = " << coeff[1]
      << " c = " << coeff[2]
      << " (";
    switch (get_type(it->first))
    {
      case DIRICHLET:
        std::cout << "Dirichlet";
        break;
      case NEUMANN:
        std::cout << "Neumann";
        break;
      default:
        std::cout << "Mixed";
        break;
    }
    std::cout << ")" << std::endl;
  }
}
