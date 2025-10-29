// $Id$

#include "GaussDOS.h"
#include "Constants.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "InitFailedException.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;


GaussDOS::GaussDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _E0(0), _sigma(0.1)
{

}

void
GaussDOS::read_database(void)
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

    vector<double> bg(3,0);
    bg[0] = db.get("Eg_G", 1e3);
    bg[1] = db.get("Eg_X", 1e3);
    bg[2] = db.get("Eg_L", 1e3);

    double bandgap = *min_element(begin(bg), end(bg));

    /*cout << "The bandgap in G is equal to:   " << bg[0] << endl;
    cout << "The bandgap in X is equal to:   " << bg[1] << endl;
    cout << "The bandgap in L is equal to:   " << bg[2] << endl;

    cout << "The min bandgap is equal to:   " << bandgap << endl;*/

    reference_energy()[0] += bandgap;

    db.set_section("conductionband");
    _N0 = db.get("N0", 1e18);
  }
  else if (get_particle() == 'h')
  {
    db.set_section("valenceband");
    reference_energy()[0] = db.get("E_v", 0.0);

    _N0 = db.get("N0", 1e18);
  }
}

void
GaussDOS::do_init(void)
{
  if (get_particle() == 'e')
  {
    if (has_parameter("level"))
    {
      get_parameter("level", reference_energy()[0]);
    }
    get_parameter("sigma", _sigma);
    get_parameter("N0", _N0);
  }
  else
  {
    if (has_parameter("level"))
    {
      get_parameter("level", reference_energy()[0]);
    }
    get_parameter("sigma", _sigma);
    get_parameter("N0", _N0);
  }
  effective_mass()[0] = 1.0;

  effective_dos() = _N0;
  total_state_density() = _N0;
}

double GaussDOS::erfc(double x) const
{
  // TODO ec is set to larger values than before. What would be 'correct'
  // values??

  if (fabs(x)<=4)
  {
    const double ec = 1e-7;
    double e;

    int n = 0;
    double Sn = 2.0*x / sqrt(M_PI);
    double S = 0.0;

    do
    {
      S += Sn;
      e = 100 * fabs(Sn / S);	

      Sn *= -1.0 * x * x * (2.0*n + 1.0) / ((n+1.0)*(2.0*n + 3.0));
      n++;
    } while (e > ec);
    return 1.0-S;
  }
  else
  {
    const double a1 = 0.254829592;
    const double a2 = -0.284496736;
    const double a3 = 1.421413741; 
    const double a4 = -1.453152027;
    const double a5 = 1.061405429;
    const double p = 0.3275911;

    if (x > 0)
    {
      double t = 1 / (1 + p*x); 	
      return ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * exp(-1*x*x);
    }
    else
    {
      x *= -1.0;
      double t = 1 / (1 + p*x);
      return 2 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * exp(-1*x*x);
    }
  }
}

double GaussDOS::inverfc(double x) const
{
  double S = 0.0;

  if ((x >= 2.0) || (x <= 0.0))
  { }
  else if (x>=0.003 && x<=1.997)
  {
    const double ec = 1e-3;
    double e;
	
    double p;
    p = 0.5 * sqrt(M_PI) * (x - 1.0);

    double f;
    f = p;

    unsigned int n = 0;
    unsigned int max_n = 100;

    long double cn;
    vector<long double> c;
    c.reserve(max_n);
    c.push_back(1.0);

    do
    {
      double Sn = c[n] * f / (2.0*n + 1.0);
      S += Sn;
      e = 100.0 * Sn / S;

      n++;
      f *= p*p;

      cn = 0.0;
      for (unsigned int m = 0; m <= n-1; m++) 
      {
        cn += (c[m]/(2.0*m+1.0))*(c[n-1-m]/(m+1.0));
      }
      c.push_back(cn);
      //cout<<"c["<<n<<"] = "<<cn<<endl;	
    } while ((e > ec) && (n < max_n));

    if (x < 1.0) 
    {
      S *= -1;
    }
  }
  else //use the asymptotic expansion
  {
    if (x < 1)
    {
      S = sqrt( - log(x) - 0.5 * log(- M_PI * log(x) - 0.5 * log(- M_PI * log(x))));
    }
    else
    {
      x = 2 - x;
      S = - sqrt( - log(x) - 0.5 * log(- M_PI * log(x) - 0.5 * log(- M_PI * log(x))));
    }
  }

  return(S);
}

inline
double GaussDOS::H(double x) const
{
  return sqrt(2.0) * inverfc(exp(-0.5 * x * x)) / x;
}

inline
double GaussDOS::K(double x, double h) const
{ 
  return 2.0*(1.0 - exp(0.5 * x * x * (1.0 - h * h)) * (h * sqrt(2.0/M_PI))/x);
}

void
GaussDOS::calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
    double kT, double, const Elem* elem, const Point& p) const
{

  double dens, der, der2;

  double N0 = get_effective_dos(elem, p);

  double ref_en = get_reference_energy()[0];
  if (get_particle() == 'h') ref_en *= -1.0;
  //cout<<"ref_energy "<<get_reference_energy()[0]<<endl;
  double z = (Ef - ref_en - Epot) / kT;

  double s = _sigma / kT;
  //cout<<"particle="<<get_particle()<<" Ef="<<Ef<<" E0="<<_E0<<" epot="<<Epot<<" z="<<z<<endl;

  double hs;
  hs = H(s);

  if (z <= -1.0*s*s) 
  {
    double espf;
    double ks;
    ks = K(s, hs);
    espf = exp( ks * (z + s*s));
    dens = N0 * exp(z + 0.5 *s*s) / (espf + 1.0);
    if (result.size() > 1)
      der = dens * (1.0 - (ks * espf / (espf + 1.0) ) ) / kT;
    if (result.size() > 2)
      der2 = der * (1.0 - (ks * espf / (espf + 1.0) ) ) / kT + dens * (ks/kT/kT) * espf / (espf+1.0) /  (espf+1.0);


    //cout<<"s="<<s<<"<z="<<z<<" dens="<<dens<<" der="<<der<<endl;	
  }
  else
  {
    dens = N0 * 0.5 * erfc(-1.0 * hs * z / (s * sqrt(2.0)) );
    if (result.size() > 1)
      der = N0 * hs * exp(-0.5 * hs * hs * z * z / (s*s)) / (s * sqrt(2.0*M_PI) * kT);
    if (result.size() > 2)
      der2 = der * (- hs * hs * z / (s*s))/kT;

  }

  result[0] = dens;
  if (result.size() > 1)
    result[1] = der;
  if (result.size() > 2)
    result[2] = der2;

}

void
GaussDOS::do_print_info(void)
{
  ostringstream os;
  os << "Gaussian center = " <<  reference_energy()[0] << " eV\n"
       << "N0 = " << _N0 << " cm^-3"
       << ", sigma = " << _sigma << " eV\n";
  Messages::info(os.str());
}

