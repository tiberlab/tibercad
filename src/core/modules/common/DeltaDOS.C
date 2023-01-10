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
DeltaDOS::calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
    double kT, double , const Elem* , const Point& ) const
{

  double dens, der, der2;

  double expf = exp(-(Ef - get_reference_energy()[0] - Epot) / kT);
  dens = _N0 / (1.0 + expf);

  if (result.size() > 1)
  {
    der = dens;
    der /= kT * (1.0 + expf);
    der *= expf;
  }
  if (result.size() > 1)
  {
    der2 = der + dens/kT * expf / (1.0 + expf);
    der2 *= expf;
    der2 /= kT * (1.0 + expf);
  }

  result[0] = dens;
  if (result.size() > 1)
    result[1] = der;
  if (result.size() > 2)
    result[2] = der2;
}
