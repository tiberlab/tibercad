
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
  const Database& db = get_database();

  db.set_section("valenceband");

  double mp = db.get("m_dos", 1.0);

  db.set_section("conductionband");

  double mn = db.get("m_dos", 1.0);

  effective_mass()[0] = mn+mp;

  db.set_section("bandgap");
  _energy = db.get("Eg_G", 1);
  _R = db.get("exciton_binding_energy", _R);
}


void
ExcitonGaussDOS::do_init(void)
{

  _R = get_option("R", _R);

  if (_R < 0)
    throw InitFailedException("Exciton binding energy R cannot be negative");

  _sigma = get_option("sigma", _sigma);
  _J = get_spin();
  _J = get_option("spin", _J);
  effective_dos() = get_option("N0", effective_dos()); // get the effective DOS if not specified from carriers

  if (has_option("N0"))
    fixed_dos() = true;

  _energy = get_option("level", _energy);
  _energy = get_option("exciton_energy", _energy);

  if (_R > _energy)
    throw InitFailedException("Exciton binding energy R cannot be greater than band gap");

  reference_energy()[0] = _energy - _R;

  
  /*
  std::ofstream of("exciton_density.dat");
  of << "# E0 = " << reference_energy()[0] << ",  sigma = " << _sigma << "\n";
  for (double E = reference_energy()[0] - 10*_sigma; E < reference_energy()[0] - 4*_sigma; E += _sigma / 100)
  {
    std::pair<double, double> res = calculate_density_and_derivative(E, 0, 0.0258, 0.0258);
    of << E << " " << res.first << " " << res.second << "\n";
  }
  */
  

}

void
ExcitonGaussDOS::do_reinit(void)
{

}
void
ExcitonGaussDOS::do_print_info(void)
{
  ostringstream os;
  os << "Exciton energy = " <<  _energy << " eV\n"
     << "Binding energy = " << _R << " eV\n";
  Messages::info(os.str());
}



void
ExcitonGaussDOS::calculate_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  return calculate_density_and_derivative(den_and_der, Ef, Epot, kT, kTlattice);
}

void
ExcitonGaussDOS::calculate_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot, double kT, double kTlattice) const
{
  double dens, der, der2;
  double N0 = (2.0*_J + 1.0) * get_effective_dos();

  std::vector<double> x, y, dy;
  x.reserve(_order + 1);
  y.reserve(_order + 1);
  dy.reserve(_order + 1);

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
    double f = _f(val, -Ef, kT);
    double gauss = exp( - val * val);
    double g = fac * gauss * f;

    y.push_back(g);
    dy.push_back(- beta * g * (f+1.0));
  }

  dens = N0 * _trapez(x, y);
  if (den_and_der.size() > 1)
    der  = -N0 * _trapez(x, dy);
  if (den_and_der.size() > 2)
    der2 = 0; //TODO

/*
  double s = 0.5 * _sigma * _sigma / kT;
  double arg = (s - energy - Ef)/kT;
  dens = (2.0*_J + 1.0) * N0 * exp(arg);
  der =  -dens / kT;
*/

  den_and_der = {dens};
  if (den_and_der.size() > 1)
    den_and_der = {dens, der};
  if (den_and_der.size() > 2)
    den_and_der = {dens, der, der2};
}


double
ExcitonGaussDOS::_f(double x, double Ef, double kT) const
{
  double energy = get_reference_energy()[0];

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














