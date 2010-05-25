// $Id$

#include "Trap.h"

//TIBER_MODULE(Trap, trap)

using namespace std;

Trap::Trap(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _density(0.0),
  _type(NEUTRAL),
  _particle('e'),
  _level(0.0)
{
  string type = get_option("type", "");
  if (type == "eNeutral")
  {
    _particle = 'e';
    set_type(NEUTRAL);
  }
  else if (type == "hNeutral")
  {
    _particle = 'h';
    set_type(NEUTRAL);
  }
  else if (type == "donor")
  {
    _particle = 'e';
    set_type(CHARGED);
  }
  else if (type == "acceptor")
  {
    _particle = 'h';
    set_type(CHARGED);
  }
  else if (type == "fixed_charge")
  {
    set_type(FIXED);
  }
}



void
Trap::do_init(void)
{
  if (_particle == 'h')
    _energy_reference = 'v';
  else
    _energy_reference = 'c';

  string tmp("c");
  tmp[0] = _energy_reference;
  tmp =  get_option("reference", tmp);
  _energy_reference = tmp[0];

  get_parameter("Nt", _density);
  get_parameter("Et", _level);
}


inline
double
Trap::_trap_level(void) const
{
  double ref;
  switch (_energy_reference)
  {
    case 'v':
      ref = _Ev + _level;
      break;

    case 'm':
      ref = 0.5 * (_Ev + _Ec) + _level;
      break;

    default:
      ref = _Ec - _level;
      break;
  }

  return ref;
}


double
Trap::get_ionized_density(void) const
{
  double dens = _density;

  if (_type != FIXED)
  {
    double f;
    double arg = _trap_level() - _fermi_level;
    double g = 1;
    double Nt = _density;

    switch (_particle)
    {
      case 'h':
        f = 1.0 / (1.0 + g * exp(-arg / _kT));
        break;

      case 'e':
      default:
        Nt = -Nt;
        f = 1.0 / (1.0 + exp(arg / _kT) / g);
        break;
    }

    dens = Nt * f;

    if (_type == CHARGED)
      dens -= Nt;
  }

  return dens;
}



double
Trap::get_ionized_density_derivative(void) const
{
  double deriv = 0.0;

  if (_type != FIXED)
  {
    double arg = _trap_level() - _fermi_level;
    double g = 1;
    double Nt = _density;
    double expfac;

    switch (_particle)
    {
      case 'h':
        expfac = g * exp(-arg / _kT);
        break;

      case 'e':
      default:
        expfac = exp(arg / _kT) / g;
        break;
    }

    double denom = 1.0 + expfac;

    deriv = -Nt / _kT * expfac / (denom * denom);
  }

  return deriv;
}
