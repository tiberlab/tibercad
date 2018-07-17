#include "ConstantDOS.h"
#include "Constants.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "InitFailedException.h"

#include "TiberModule.h"

using namespace std;

ConstantDOS::ConstantDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _Ewidth(0.0), _N0(1.0)
{

}

void
ConstantDOS::read_database(void)
{
  // when reading from the database, we use the same data
  // as for kp
  const Database& db = get_database();

  if (get_particle() == 'e')
  {
    // TODO: should bowing be applied to Eg(T) or Eg(0) ?

    db.set_section("valenceband");
    reference_energy()[0] = db.get("E_v", 0.0);

    db.set_section("bandgap");
    double bandgap = db.get("Eg_G", 1e3);

    reference_energy()[0] += bandgap;
  }
  else if (get_particle() == 'h')
  {
    db.set_section("valenceband");
    reference_energy()[0] = db.get("E_v", 0.0);
  }
}

void
ConstantDOS::do_init(void)
{
  if (has_parameter("level"))
    get_parameter("level", reference_energy());

  get_parameter("Ewidth", _Ewidth);
  get_parameter("N0", _N0);

  effective_mass()[0] = 1.0;

  effective_dos() = _N0;
  total_state_density() = _N0;
}

void
ConstantDOS::calculate_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  return calculate_density_and_derivative(den_and_der, Ef, Epot, kT, kTlattice);
}

void
ConstantDOS::calculate_density_and_derivative(std::vector<double>& den_and_der, double Ef, double Epot, double kT, double kTlattice) const
{
  double dens, der, der2, Emin, Emax, exp_num, exp_den, ref_en;
  //cout<<"ref_energy "<<reference_energy()<<endl;

  ref_en = get_reference_energy()[0];
  if (get_particle() == 'h') ref_en *= -1.0;

  if (get_particle() == 'e')
  {
    Emin = ref_en;
    Emax = Emin + _Ewidth; 
  }
  else
  {
    Emax = ref_en;
    Emin = Emax + _Ewidth;
  }
  
  double zmin = (Ef - Emin - Epot) / kT;
  double zmax = (Ef - Emax - Epot) / kT;
  
  if (get_particle() == 'e')
  {
    exp_num = exp(zmin);
    exp_den = exp(zmax);
  }
  else
  {
    exp_num = exp(zmax);
    exp_den = exp(zmin);
  }

  double C = _N0 / _Ewidth;

  dens = kT * C * log((exp_num + 1)/(exp_den +1));
  if (den_and_der.size() > 1)
    der = C * ( exp_num / (exp_num + 1) - exp_den / (exp_den + 1));
  if (den_and_der.size() > 2)
    der2 = C/kT * (std::pow(exp_num/(exp_num + 1),2) - std::pow(exp_den/(exp_den + 1),2));

  
  den_and_der[0] = dens;
  if (den_and_der.size() > 1)
    den_and_der[1] = der;
  if (den_and_der.size() > 2)
    den_and_der[2] = der2;

}
