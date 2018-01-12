// $Id: KPBulkDOS.C 3347 2012-06-15 15:00:26Z maufder $

#include "KPBulkDOS.h"
#include "Constants.h"
#include "Material.h"
#include "TiberMath.h"
#include "Database.h"
#include "InitFailedException.h"
#include "DDsemiconductor.h"
#include "SimulationOptions.h"
#include "SimulationInterface.h"
#include "Messages.h"

#include "TiberModule.h"

#include "elem.h"


using namespace std;

KPBulkDOS::KPBulkDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _recompute(false)
{
}




void
KPBulkDOS::prepare_submodels(void)
{
  ModelOptions opts(get_options());
  opts.delete_all_submodels();
  opts.delete_option("strain");
  opts.delete_option("strain_simulation");
  opts.delete_option("particle");
  _bulk_model = DDsemiconductor::create(get_material(), opts);

  if (_bulk_model == NULL)
    throw InitFailedException("Cannot create KdotP bulk model for material "
        + get_material()->get_name());

  add_submodel("kp", _bulk_model);
}



void
KPBulkDOS::do_init(void)
{

  /*
  _ref_energies.resize(0);
  get_parameter("level", reference_energy());
  _ref_energies.push_back(reference_energy());

  get_option("levels", _ref_energies);

  get_option("dos_mass", _dos_mass);


  get_option("degeneracy", _degeneracy);


  if (_dos_mass.size() != _ref_energies.size())
    throw InitFailedException("Bulk DOS model for material " + get_material()->get_name() +
        ": number of energy levels different from number of DOS masses");

  if (_degeneracy.size() != _ref_energies.size())
    throw InitFailedException("Bulk DOS model for material " + get_material()->get_name() +
        ": number of energy levels different from number of degeneracies");
  */

  _dos_factor = pow(2.0 * M_PI *
        Constants::me / (Constants::h * Constants::h) *
        Constants::e, 1.5) / 1e6;

  string strain_simul = get_option("strain_simulation", "");
  strain_simul = get_option("strain", strain_simul);
  _strain_if.set_simulation(strain_simul);

  // TODO temperature interface !!

  _bulk_model->set_strain(Tensor2Sym(0.0));
  _bulk_model->set_temperature(SimulationOptions::T);
  _solve_kp();

  if (_strain_if.has_simulation())
    _recompute = true;

}



void
KPBulkDOS::do_reinit(const Elem* elem)
{
  if (!_recompute)
    return;

  Cache::Data data;
  if (_cache.get_data(elem,
      _strain_if.get_simulation()->get_solve_sequence_number(),
      data))
  {
    effective_dos() = data.eff_dos;
    band_edge() = data.ref_energy;
    reference_energy() = data.ref_energy;
    effective_mass() = data.dos_mass;
    dos_mass() = data.dos_mass;
  }
  else
  {
    Tensor2Sym strain(0);
    _strain_if.get_crystal_strain(elem, elem->centroid(), strain);
    _bulk_model->set_strain(strain);

    // it wants temperature in K
    _bulk_model->set_temperature(SimulationOptions::T);

    _solve_kp();

    data.strain_seq_num = _strain_if.get_simulation()->get_solve_sequence_number();
    data.eff_dos = effective_dos();
    data.ref_energy = reference_energy();
    data.dos_mass = effective_mass();
    _cache.set_data(elem, data);
  }
}


void
KPBulkDOS::do_print_info(void)
{
  ostringstream os;
  os << "band edge :";
  for (int i = 0; i < band_edge().size(); ++i)
    os << "  " << band_edge()[i];
  Messages::info(os.str());

  os.str("");
  os << "eff. mass :";
  for (int i = 0; i < band_edge().size(); ++i)
    os << "  " << dos_mass()[i];
  Messages::info(os.str());

  os.str("");
  os << "degeneracy:";
  for (int i = 0; i < band_edge().size(); ++i)
    os << "  " << degeneracy()[i];
  Messages::info(os.str());
}


void
KPBulkDOS::_solve_kp(void)
{
  double kT = SimulationOptions::T * Constants::k_B;

  // calculate conduction and valence band data
  if (get_particle() == 'e')
  {
    _bulk_model->calculate_conduction_band_extremum();

    const std::vector<DDsemiconductor::band_extremum>& cbs =
        _bulk_model->get_conduction_band_energy_mass();

    // get minimum

    band_edge().resize(cbs.size());
    dos_mass().resize(cbs.size());
    degeneracy().resize(cbs.size());
    band_edge()[0] = cbs[0].energy;
    dos_mass()[0] = cbs[0].mass_DOS;
    degeneracy()[0] = cbs[0].degeneracy;

    for (unsigned int i = 1; i < cbs.size(); i++)
    {
      band_edge()[i] = cbs[i].energy;
      dos_mass()[i] = cbs[i].mass_DOS;
      degeneracy()[i] = cbs[i].degeneracy;
      //if (cbs[i].energy < band_edge()[0])
      //{
      //  swap(band_edge()[i], band_edge()[0]);
      //  swap(dos_mass()[i], dos_mass()[0]);
      //  swap(degeneracy()[i], degeneracy()[0]);
      //}
    }
  }
  else if (get_particle() == 'h')
  {
    _bulk_model->calculate_valence_band_extremum();

    const std::vector<DDsemiconductor::band_extremum>& vbs =
      _bulk_model->get_valence_band_energy_mass();

    // get maximum

    band_edge().resize(vbs.size());
    dos_mass().resize(vbs.size());
    degeneracy().resize(vbs.size());
    band_edge()[0] = vbs[0].energy;
    dos_mass()[0] = vbs[0].mass_DOS;
    degeneracy()[0] = vbs[0].degeneracy;

    for (unsigned int i = 1; i < vbs.size(); i++)
    {
      band_edge()[i] = vbs[i].energy;
      dos_mass()[i] = vbs[i].mass_DOS;
      degeneracy()[i] = vbs[i].degeneracy;

      //if (vbs[i].energy > band_edge()[0])
      //{
      //  swap(band_edge()[i], band_edge()[0]);
      //  swap(dos_mass()[i], dos_mass()[0]);
      //  swap(degeneracy()[i], degeneracy()[0]);
      //}
    }

    /*
    double delta_max = 4.0 * kT;
    double tmp = 0;
    // include other bands
    for (unsigned int i = 0; i < vbs.size(); i++)
    {
      double delta = band_edge()[0] - vbs[i].energy;
      if (delta < delta_max)
        tmp += vbs[i].degeneracy * std::pow(vbs[i].mass_DOS, 1.5)
          * std::exp(-delta / kT);
    }
    dos_mass()[0] = std::pow(tmp / degeneracy()[0], 2.0 / 3.0);
    */
  }

  effective_dos() = _dos_factor * _degeneracy[0] * std::pow(kT * _dos_mass[0], 1.5);
  reference_energy() = band_edge();
  effective_mass() = dos_mass();
}





std::pair<double, double>
KPBulkDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  double density = 0.0;
  double derivative = 0.0;

  const double arg_max = 150;
  const double arg_min = -50;
  const double min_dens = 1e-64;

  _th_el_power = 0;

  for (int i = 0; i < _ref_energies.size(); ++i)
  {
    double dens, der;

    double energy = _ref_energies[i];

    if (get_particle() == 'e')
    {
      //double T = kTlattice / Constants::k_B;
      //energy += _gap_params[i].Eg0 -
      //    _gap_params[i].varshni_a * T * T / (T + _gap_params[i].varshni_b);
    }
    else if (get_particle() == 'h')
      energy = -energy;

    double arg = (Ef - energy - Epot) / kT;

    double Neff = _dos_factor * _degeneracy[i] * std::pow(kT * _dos_mass[i], 1.5);

    if (arg < arg_min)
    {
      dens = exp(arg);
      der = dens;
    }
    else if (arg < arg_max)
    {
      dens = TiberMath::fermidirac_half(arg);
      der = TiberMath::fermidirac_mhalf(arg);
    }
    else
    {
      dens = 2.0 * M_2_SQRTPI / 3.0 * std::pow(arg, 1.5);
      der = M_2_SQRTPI * std::sqrt(arg);
    }

    dens *= Neff;

    // is this needed???
    if (dens > min_dens)
      der *= Neff / kT;
    else
      der = 0.0;

    density += dens;
    derivative += der;


    // calculate thermoelectric power
    double temp = kT / Constants::k_B;
    double Pth = dens / der * 1.5 / temp + Constants::k_Boltzmann * (1 - arg);

    _th_el_power += dens * Pth;

  }

  _th_el_power /= density;


  return make_pair(density, derivative);
}

std::pair<double, double>
KPBulkDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice) const
{
  return make_pair(0, 0);
}

