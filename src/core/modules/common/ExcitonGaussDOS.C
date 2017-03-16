
#include "ExcitonGaussDOS.h"
#include "Constants.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "InitFailedException.h"
#include "SimulationOptions.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;

ExcitonGaussDOS::ExcitonGaussDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _order(64),
  _R(0.2),
  _J(0.0),
  _sigma(0.2)
{
}



void
ExcitonGaussDOS::read_database(void)
{
  // when reading from the database, we use the same data
  // as for kp
  const Database& db = get_database();

  db.set_section("valenceband");

  double mp = db.get("m_dos", 1.0);

  db.set_section("conductionband");

  double mn = db.get("m_dos", 1.0);

  effective_mass()[0] = mn+mp;
}


void
ExcitonGaussDOS::do_init(void)
{

  _R = get_option("R", _R);
  _sigma = get_option("sigma", _sigma);
  _J = get_spin();
  _J = get_option("spin", _J);
  effective_dos() = get_option("N0", effective_dos()); // get the effective DOS if not specified from carriers

  if (has_option("N0"))
    fixed_dos() = true;

  double level = get_option("level", 1.0);
  double energy = get_option("exciton_energy", 1.0);
  reference_energy()[0] = level - _R;

  if (reference_energy().size() == 3)
    reference_energy()[0] = max(reference_energy()[1], reference_energy()[2]) - _R;
  else
  {
    reference_energy().resize(2);
    reference_energy()[1] = energy;  //store exciton energy if not specified from carriers
  }

  double egap = fabs(reference_energy()[1] - reference_energy()[2]);

  if (_R < 0)
    throw InitFailedException("Exciton binding energy R cannot be negative");

  if (_R > egap)
    throw InitFailedException("Exciton binding energy R cannot be greater than band gap");

}

void
ExcitonGaussDOS::do_reinit(void)
{

}
void
ExcitonGaussDOS::do_print_info(void)
{
  ostringstream os;
  double energy = (reference_energy().size() == 3) ? fabs(reference_energy()[1]-reference_energy()[2]) - _R : reference_energy()[1];
  os << "Exciton energy = " <<  energy << " eV\n"
     << "Binding energy = " << _R << " eV\n";
  Messages::info(os.str());
}



std::pair<double, double>
ExcitonGaussDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  return calculate_density_and_derivative(Ef, Epot, kT, kTlattice);
}

std::pair<double, double>
ExcitonGaussDOS::calculate_density_and_derivative(double Ef, double Epot, double kT, double kTlattice) const
{
  double dens, der;
  double N0 = (2.0*_J + 1.0) * get_effective_dos();

  std::vector<double> x, y, dy;

  double s = 1.0 / sqrt(2.0);
  double xM = 4.0 * s;
  double x0 = -xM;
  double step = (xM - x0) / _order;

  for (unsigned int i=0; i <= _order; i++)
    x.push_back(i*step + x0);

  double fac = 1.0 / sqrt(M_PI);
  double beta = 1.0 / kT;

  for (auto val : x)
  {
    double f = _f(val, Ef, kT);
    double gauss = exp( - val * val);
    double g = fac * gauss * f;

    y.push_back(g);
    dy.push_back(- beta * g * (f+1.0));
  }

  dens = N0 * _trapez(x, y);
  der  = N0 * _trapez(x, dy);

/*
  double s = 0.5 * _sigma * _sigma / kT;
  double arg = (s - energy - Ef)/kT;
  dens = (2.0*_J + 1.0) * N0 * exp(arg);
  der =  -dens / kT;
*/

  return make_pair(dens, der);
}


double
ExcitonGaussDOS::_f(double x, double Ef, double kT) const
{
  const vector<double>& refenergy = get_reference_energy();
  double energy = (refenergy.size() == 3) ? fabs(refenergy[1]-refenergy[2]) - _R : refenergy[1];

  double x0 = energy / _sigma / sqrt(2.0);
  double xf = Ef / _sigma / sqrt(2.0);

  if (x <= -(x0+xf))
    return 0.0;

  double beta = 1.0/kT;
  double arg = beta * _sigma * sqrt(2.0) * (x + x0 + xf);
  double exponential = exp(arg);
  double f = 1.0 / (exponential - 1.0);

  return f;
}

double
ExcitonGaussDOS::_trapez(std::vector<double>& x, std::vector<double>& y) const
{
  assert(x.size() == y.size());

  double integral = 0.0;

  for (unsigned int i = 1; i<x.size(); i++)
    integral += (y[i-1]+y[i])*(x[i]-x[i-1]) * 0.5;

  return integral;
}














