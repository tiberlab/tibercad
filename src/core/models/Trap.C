// $Id$

#include "Trap.h"
#include "Particle.h"
#include "DensityOfStates.h"
#include "TiberMath.h"

#include "fstream"

using namespace std;

bool Trap::_coupled;



Trap::Trap(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _density(0.0),
  _profile(nullptr),
  _type(NEUTRAL),
  _particle('e'),
  _level(0.0),
  _energy_reference('m'),
  _sigma_n(1e-15),
  _sigma_p(1e-15),
  _tau_n(0.0),
  _tau_p(0.0),
  _e_vth(1e7),
  _h_vth(1e7),
  _gen_TC(0.0),
  _gen_VT(0.0),
  _dos(nullptr),
  _ext_dens_sim(nullptr)
{
  string type = get_option("type", "");
  if (type == "eNeutral")
  {
    _particle = 'e';
    set_type(NEUTRAL);
  }
  else if (type == "hNeutral")
  {
    _particle = 'h';
    set_type(NEUTRAL);
  }
  else if (type == "donor")
  {
    _particle = 'e';
    set_type(CHARGED);
  }
  else if (type == "acceptor")
  {
    _particle = 'h';
    set_type(CHARGED);
  }
  else if (type == "fixed_charge")
  {
    set_type(FIXED);
  }
  else
  {
    throw InitFailedException("Unknown trap type '" + type + "'");
  }
}


Trap::~Trap(void)
{
  destroy(_dos);
  delete _profile;
}


void
Trap::prepare_submodels(void)
{
  if (get_options().has_submodel("density_of_states"))
  {
    ModelOptions::submodel_iterator it(get_options().submodels_begin("density_of_states"));
    _dos = DensityOfStates::create(it->second);
    add_submodel("dos", _dos);
  }

}


void
Trap::do_init(void)
{
  string tmp("m");
  tmp =  get_option("reference", tmp);
  _energy_reference = tmp[0];

  if (get_options().has_submodel("profile"))
  {
    _profile = ExternalProfile::create(
        get_options().submodels_begin("profile")->second);
    _density = 1.0;
  }

  get_parameter("Nt", _density);
  get_parameter("Et", _level);

  get_parameter("sigma_n", _sigma_n);
  get_parameter("sigma_p", _sigma_p);

  get_parameter("tau_n", _tau_n);
  get_parameter("tau_p", _tau_p);

  // TODO to be obtained from DOS model (?)
  get_parameter("thermal_velocity_n", _e_vth);
  get_parameter("thermal_velocity_p", _h_vth);

  get_parameter("trap_to_cb_rate", _gen_TC);
  get_parameter("vb_to_trap_rate", _gen_VT);

  /* for testing
  double Ec = 1,  Ev = 0;
  string name = "trap_";
  name += _particle;
  name += ".dat";
  ofstream of(name.c_str());
  double h = 1.0 / 10000;
  of << "# " << Ec << " " << Ev << "\n";
  for (unsigned int i = 0; i < 10000; i++)
  {
    double phi = i * h;
    double Efn = 0;
    double Efp = 0;

    double n = 1e19*exp(-(Ec-phi) / 0.026);
    double p = 1e19*exp((Ev-phi) / 0.026);
    Particle el(-1, n, Efn, 0.026);
    Particle hl(1, p, Efp, 0.026);
    vector<double> der;
    set_energies(Ec-phi, Ev-phi);
    double f = get_ionized_density_and_derivative(el, hl, der);
    of << phi << " " << n << " " << p  << " " <<
        f << " " << der[0] << " " << der[1] << " " << der[2] << " " << der[3] <<  "\n";
  }
  of.close();
  */
  _coupled = false;
  if (_type == FIXED ) 
  {
    string sim = get_option("ext_dens_simulation", "");

    if ( sim != "" )
      _ext_dens_sim = SimulationInterface::find_simulation(sim);

    
    if ( ( _ext_dens_sim == nullptr ) && ( sim != "" ) )
    {
      std::string msg("External density simulation " + std::string(sim) + " not found");
      throw InitFailedException(msg);
    }

    if (_ext_dens_sim != nullptr)
    {
      _eDensity = _ext_dens_sim->get_solution_id("eDensity");
      _hDensity = _ext_dens_sim->get_solution_id("hDensity");

      _coupled = true; 
    }
  }
}


inline
double
Trap::_trap_level(void) const
{
  double ref;
  switch (_energy_reference)
  {
    case 'v':
      ref = _Ev + _level;
      break;

    case 'm':
      ref = 0.5 * (_Ev + _Ec) + _level;
      break;

    default:
      ref = _Ec - _level;
      break;
  }

  return ref;
}



double
Trap::get_ionized_density_and_derivative(const Elem* elem, const Point& p,
    const Particle& el, const Particle& hl, std::vector<double>& derivatives) const
{
  double dens = _density;
  if (_profile != nullptr)
    dens *= _profile->get_data(elem, p);

  double Nt = dens;

  derivatives.resize(5);
  derivatives[0] = derivatives[1] = derivatives[2] = derivatives[3] = derivatives[4] = 0.0;

  if ((_type != FIXED) && (dens > 0.0))
  {
    double f_e, f_h;
    double deriv_e, deriv_h;
    double g = 1;

    double kT_e = el.kT();
    double kT_h = hl.kT();
    double arg_e = _trap_level() + el.fermi_level();
    double arg_h = _trap_level() + hl.fermi_level();

    double Cn = 1.0 /_tau_n / Nt;
    double Cp = 1.0 /_tau_p / Nt;

    //cout << Cn << Cp << endl;

    if (_tau_n == 0.0)
    {
      Cn = _sigma_n * _e_vth;
    }

    if (_tau_p == 0.0)
    {
      Cp = _sigma_p * _h_vth;
    }


    double n = el.density();
    double p = hl.density();

    // occupations in terms of electrons
    if (_dos == nullptr)
    {
      std::pair<double, double> occ_e(Distributions::fermi_dirac(-arg_e, kT_e));
      f_e = occ_e.first;
      deriv_e = occ_e.second;

      std::pair<double, double> occ_h(Distributions::fermi_dirac(-arg_h, kT_h));
      f_h = occ_h.first;
      deriv_h = occ_h.second;


      double f;
      if (f_e < 1e-12)
      {
        f_e = 1e-12;
        deriv_e = 0;
      }
      if (f_h > (1.0 - 1e-12))
      {
        f_h = 1.0 - 1e-12;
        deriv_h = 0;
      }

      switch (_particle)
      {
        case 'h':
        {
          double gc = (1.0 - f_e) / f_e;
          double gv = f_h / (1.0 - f_h);

          //double nom = Cp * p + Cn * n * gc - (_gen_VT - _gen_TC) / Nt;
          //double denom = Cn * n * (1 + gc) + Cp * p * (1 + gv);
          double nom = Cp * p + Cn * n * gc + _gen_TC / Nt;
          double denom = Cn * n * (1 + gc) + Cp * p * (1 + gv) + (_gen_VT + _gen_TC) / Nt;

          // refactorized to prevent numerical problems
          //double nom = ((Cp * p + _gen_TC / Nt) * f_e + Cn * n * (1 - f_e)) * (1 - f_h);
          //double denom = (Cn * n  + f_e * (_gen_VT + _gen_TC) / Nt) * (1 - f_h)
          //    + Cp * p * f_e;

          f = nom / denom;

          // in some cases f_e ~ 0 or f_h ~ 1, and both nom and denom may become 0
          if (nom == 0.0)
            f = 0.0;
          else if (denom == 0.0)
          {
            f = 1.0;
            denom = 1e-12;
          }

          double dfdp = Cp * (1.0 - f * (1 + gv)) / denom;
          double dfdn = Cn * (gc - f * (1 + gc)) / denom;

          double dfdEfn = -(deriv_e / f_e) * ((1 - f) / f_e) * Cn * n / denom;
          double dfdEfp = -(f / (1 - f_h)) * (deriv_h / (1 - f_h)) * Cp * p / denom;

          // this is somewhat crude ...
          if (denom == 0.0)
            dfdp = dfdn = dfdEfn = dfdEfp = 0.0;

          derivatives[0] = Nt * dfdn;
          derivatives[1] = Nt * dfdp;
          derivatives[2] = Nt * dfdEfn;
          derivatives[3] = Nt * dfdEfp;

          break;
        }

        case 'e':
        default:
        {
          double gc = (1.0 - f_e) / f_e;
          double gv = f_h / (1.0 - f_h);

          //double nom = Cn * n + Cp * p * gv + (_gen_VT - _gen_TC) / Nt;
          //double denom = Cn * n * (1 + gc) + Cp * p * (1 + gv);
          double nom = Cn * n + Cp * p * gv + _gen_VT / Nt;
          double denom = Cn * n * (1 + gc) + Cp * p * (1 + gv) + (_gen_VT + _gen_TC) / Nt;

          f = nom / denom;

          // in some cases f_e ~ 0 or f_h ~ 1, and both nom and denom may become 0
          if (nom == 0.0)
            f = 0.0;
          else if (denom == 0.0)
          {
            f = 1.0;
            denom = 1e-12;
          }

          double dfdn = Cn * (1.0 - f * (1 + gc)) / denom;
          double dfdp = Cp * (gv - f * (1 + gv)) / denom;

          double dfdEfn = (deriv_e / f_e) * (f / f_e) * Cn * n / denom;
          double dfdEfp = ((1 - f) / (1 - f_h)) * (deriv_h / (1 - f_h)) * Cp * p / denom;

          // this is somewhat crude ...
          if (denom == 0.0)
            dfdp = dfdn = dfdEfn = dfdEfp = 0.0;

          Nt = -Nt;

          derivatives[0] = Nt * dfdn;
          derivatives[1] = Nt * dfdp;
          derivatives[2] = Nt * dfdEfn;
          derivatives[3] = Nt * dfdEfp;

          break;
        }
      }

      //cout << "f = " << f << endl;
      dens = Nt * f;
    }  // end _dos == nullptr

    else // _dos != nullptr
    {
      double f, deriv;
      double level = _trap_level();
      //double Nt = _density;

      switch (_particle)
      {
        case 'h':
        {
          // set the reference energy
          // NOTE: need to take away electrostatic energy because it is
          // added internally in the DOS model
          _dos->set_reference_energy(- level - _phi);
          //cout<<"- hl.fermi_level() = " << - hl.fermi_level() << " _phi = " << _phi << endl;
          //cout<<"hl.fermi_level() = " << hl.fermi_level() << endl ;
          std::pair<double, double> result(
              _dos->get_occupied_density_and_derivative(hl.fermi_level(), _phi, kT_h));
          f = result.first;
          deriv = -Nt * result.second;
          //cout<<"f_h = " << f << " deriv_h = " << deriv << endl;
          break;
        }
        case 'e':
        default:
        {
          Nt = -Nt;
          // it needs the fermi level shifted by trap_level
          _dos->set_reference_energy(level + _phi);
          //cout<<"level + _phi = " << level + _phi << endl;
          
          std::pair<double, double> result(
              _dos->get_occupied_density_and_derivative(-el.fermi_level(), - _phi, kT_e));
          f = result.first;
          deriv = Nt * result.second;
          //cout<<"f_e = " << f << " deriv_e = " << deriv << endl;
          break;
        }
      }

      derivatives[4] = deriv;
      dens = Nt * f;
    } // end _dos != nullptr

    if (_type == CHARGED)
      dens -= Nt;
  }

  
  if ( ( _type == FIXED ) && ( _ext_dens_sim != nullptr ) )
  {
    
    if ( _ext_dens_sim->is_solved() && _coupled )
    {
      double extn, extp = 0.0;

      _coupled = false;

      _ext_dens_sim->get_solution(elem, _eDensity, extn, p);
      _ext_dens_sim->get_solution(elem, _hDensity, extp, p);

      _coupled = true;

      dens += extp - extn;
    }
    
  }
  
  return dens;
}


