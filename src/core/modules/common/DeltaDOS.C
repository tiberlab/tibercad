// $Id$

#include "DeltaDOS.h"

#include "TiberModule.h"

using namespace std;

DeltaDOS::DeltaDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _N0(1.0)
{
}


void
DeltaDOS::do_init(void)
{
  get_parameter("level", reference_energy());
  _N0 = get_option("N0", _N0);
  effective_mass()[0] = 1.0;
  effective_dos() = _N0;
  total_state_density() = _N0;
}

std::pair<double, double>
DeltaDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  return calculate_density_and_derivative(Ef, Epot, kT, kTlattice);
}

std::pair<double, double>
DeltaDOS::calculate_density_and_derivative(double Ef, double Epot, double kT, double kTlattice) const
{
  double expf = exp(-(Ef - get_reference_energy()[0] - Epot) / kT);
  double dens = _N0 / (1.0 + expf);
  double der = dens;

  der /= kT * (1.0 + expf);
  der *= expf;

  return make_pair(dens, der);
}
