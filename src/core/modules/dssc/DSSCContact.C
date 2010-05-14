// $Id$

#include "DSSCContact.h"
#include "SimulationOptions.h"
#include "Constants.h"

using namespace std;

bool
DSSCContact::_open_circuit = true;



DSSCContact*
DSSCContact::create(const std::string& name,
    const ModelOptions& options)
{
  DSSCContact* ct = NULL;

  ct = new DSSCContact(options);

  if (ct != NULL)
    ct->set_options(options);


  if (name == "Pt")
    ct->is_cathode() = true;

  return ct;
}



void
DSSCContact::do_init(void)
{
  get_parameter("load", _res);
  get_parameter("bias", _bias);
  get_parameter("j0", _j0);
  //get_parameter("beta", _beta);
}



void
DSSCContact::calculate_current(double I, double I3)
{
  //double kT = Constants::k_B * SimulationOptions::T;
  //double upt = 0.0;
  //double A = sqrt(I3 * _Ioc / (I * _I3oc)) * exp((1 - _beta) * upt / kT);
  //double B = I / _Ioc * exp(-_beta * upt / kT);

  //_current = _j0 * (A - B);
  //_current = get_potential();
}
