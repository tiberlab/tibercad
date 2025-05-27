// $Id$

#include "BCPressure.h"
#include "SimulationOptions.h"

#include <cmath>

#include "TiberModule.h"

using namespace libMesh;


void
BCPressure::do_init(void)
{
  _relative_humidity = get_option("relative_humidity", _relative_humidity);
}


void
BCPressure::calculate(const Elem* elem, unsigned int side,
    const Point& point)
{
  double temp = SimulationOptions::temperature;
  double pressure = get_saturation_pressure(temp);
  pressure *= _relative_humidity / 100;

  set_coefficients(1, 0, pressure);
}

double
BCPressure::get_saturation_pressure(double T) const
{
    const double Tb = 373.16;     // K, boiling point of water
    const double Pb = 1013.246;   // hPa, standard atmospheric pressure

    double term1 = -7.90298 * (Tb / T - 1.0);
    double term2 = 5.02808 * std::log10(Tb / T);
    double term3 = -1.3816e-7 * (std::pow(10.0, 11.344 * (1.0 - T / Tb)) - 1.0);
    double term4 = 8.1328e-3 * (std::pow(10.0, -3.49149 * (Tb / T - 1.0)) - 1.0);
    double log10_es = term1 + term2 + term3 + term4 + std::log10(Pb);

    return 100 * std::pow(10.0, log10_es);  // e_s in hPa
}
