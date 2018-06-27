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

void
DeltaDOS::calculate_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  return calculate_density_and_derivative(den_and_der, Ef, Epot, kT, kTlattice);
}

void
DeltaDOS::calculate_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot, double kT, double kTlattice) const
{
  double dens, der, der2;

  double expf = exp(-(Ef - get_reference_energy()[0] - Epot) / kT);
  dens = _N0 / (1.0 + expf);

  if (den_and_der.size() > 1)
  {
    der = dens;
    der /= kT * (1.0 + expf);
    der *= expf;
  }
  if (den_and_der.size() > 1)
  {
    der2 = 0; //TODO
    der2 *=0; //TODO
    der2 /=0; //TODO
  }

  den_and_der.push_back(dens);
  if (den_and_der.size() > 1)
    den_and_der.push_back(der);
  if (den_and_der.size() > 2)
    den_and_der.push_back(der2);
}
