// $Id$

#include "ExponentialDOS.h"

using namespace std;

ExponentialDOS::ExponentialDOS(const ModelOptions& options) :
  DensityOfStates(options)
{
}


void
ExponentialDOS::do_init(void)
{
  get_parameter("alpha", _alpha);
}

inline
double
ExponentialDOS::get_value(double e, double E, double kT) const
{
  double num = exp(-e / _alpha);
  double den = 1 + exp((e - E) / kT);
  return num / den;
}

double
ExponentialDOS::get_occupied_density(double E, double kT) const
{
  double emin = -10 * _alpha;
  int steps = 20;
  double h = -emin / steps;
  double sum = 0;
  for (int i = 0; i <= steps; i++)
  {
    double e = emin + h * i;
    sum += get_value(e, E, kT);
  }

}


double
ExponentialDOS::get_occupied_density_derivative(double E, double kT) const
{

}
