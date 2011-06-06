// $Id$

#include "ExponentialDOS.h"

using namespace std;

ExponentialDOS::ExponentialDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _alpha(0.1)
{
}


void
ExponentialDOS::do_init(void)
{
  get_parameter("alpha", _alpha);
}

inline
double
ExponentialDOS::_get_value(double e, double E, double kT) const
{
  double num = exp(e / _alpha);
  double den = _alpha * (1 + exp((e - E) / kT));
  return num / den;
}


/*
inline
double _get_value_derivative(double e, double E, double kT) const
{
  double num = exp(e / _alpha);

}
*/


double
ExponentialDOS::get_occupied_density(double E, double kT) const
{
  double sum = 0;

  // We divide integration into 3 parts
  double g = min(5 * kT, 5 * _alpha);
  if (E - g < 0)
  {
    sum += exp((E - g) / _alpha);
    double emax = 0;

    if (E + g < 0)
    {
      if (kT != _alpha)
      {
        double tmp = (E + g) * (kT - _alpha) / (kT * _alpha);
        sum += kT / (kT - _alpha) * (1 - exp(tmp)) * exp(E/kT);
      }
      emax = E + g;
    }

    double emin = E - g;
    int steps = 25;
    double h = (emax - emin) / steps;
    for (int i = 0; i < steps; i++)
    {
      double e = emin + h * i;
      sum += h * _get_value(e, E, kT);
    }
  }
  else
    sum = 1;

  return sum;
}


double
ExponentialDOS::get_occupied_density_derivative(double E, double kT) const
{
  double der = get_occupied_density(E, kT);
  der -= 1.0 / (1 + exp(-E / kT));
  return der / _alpha;
}
