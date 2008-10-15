// $Id$


#include "FowlerNordheim.h"

#include <cmath>
#include <iostream>

using namespace std;

double
FowlerNordheim::get_emission_current(double F)
{
  double J = 0.0;
  if (fabs(F) > 1e-3)
  {
    double t_square = 1.1164;
    double nu_0 = 0.93685;
    double sqrt_wf = sqrt(_workfunction);
    double A = 1.5414e-6 / (t_square * _workfunction) * exp(9.83596 / sqrt_wf);
    double B = nu_0 * 6.8309e7 * sqrt_wf * sqrt_wf * sqrt_wf;

    // the field is opposite to that assumed in the original formula
    // due to the projection onto the normal
    J = A * F * F * exp(B / F);
  }

  return J;
}


