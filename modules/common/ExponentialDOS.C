// $Id$

#include "ExponentialDOS.h"

#include "tibercad/module/TiberModule.h"


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
  effective_mass()[0] = 1.0;
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


void
ExponentialDOS::calculate_density_and_derivative(vector<double>& result, double E, double Epot,
    double kT, double , const Elem* , const Point& ) const
{

  E -= get_reference_energy()[0] + Epot;

  double sum = 0, der = 0, der2 = 0;


  if (_alpha < 0.5 * kT)
  {
    // it's almost a delta
    double exp_E_kT = exp(-E / kT);
    sum = 1.0 / (1.0 + exp_E_kT);

    if (result.size() > 1)
      der = sum / (1 + exp_E_kT);
    if (result.size() > 2)
      der2 = - sum / (1 + exp_E_kT) + 2* sum / (1 + exp_E_kT) / (1 + exp_E_kT);

    der *= exp_E_kT / kT;
    der2 *= exp_E_kT / kT /kT;

  }
  else
  {
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

      double expm = exp(-E / kT);

      if (result.size() > 1)
      {
        der = sum - 1.0 / (1 + expm);
        der /= _alpha;
      }
      if (result.size() > 2)
      {
        der2 = sum - _alpha*(expm/(1 + expm))/kT;
        der2 /= _alpha * _alpha;
      }

    }
    else
    {
      sum = 1;
      der = 0;
      der2 = 0;
    }
  }

  result[0] = sum;
  if (result.size() > 1)
    result[1] = der;
  if (result.size() > 2)
    result[2] = der2;

}

