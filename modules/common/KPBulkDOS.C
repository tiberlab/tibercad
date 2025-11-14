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
#include "TensorOperators.h"
#include "Messages.h"

#include "libmesh/elem.h"
#include "libmesh/tensor_value.h"

#include "TiberModule.h"


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
  _strain_if = SimulationInterface::find_solution_provider(strain_simul, "Strain");

  string thermal_simul = get_option("thermal_simulation", "");
  _thermal_if = SimulationInterface::find_solution_provider(thermal_simul, "T");

  _bulk_model->set_strain(Tensor2Gen(0.0));
  _bulk_model->set_temperature(SimulationOptions::T);
  _solve_kp();

  if (_strain_if.is_valid() || _thermal_if.is_valid())
    _recompute = true;

}



void
KPBulkDOS::do_reinit(const Elem* elem)
{
  if (!_recompute)
    return;

  Cache::Data data;
  if (_cache.get_data(elem,
      _strain_if.simulation()->get_solve_sequence_number(),
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
    Tensor2Gen strain(0);

    vector<Point> p(1);
    p[0] = elem->vertex_average();
    vector<double> values(6);

    if (_strain_if.simulation()->get_solution(elem,
        _strain_if.id(), values, p))
    {
      strain(1,1) = values[0];
      strain(2,2) = values[1];
      strain(3,3) = values[2];
      strain(2,1) = strain(1,2) = values[3];
      strain(3,2) = strain(2,3) = values[4];
      strain(3,1) = strain(1,3) = values[5];
    }

    const Material* mat = get_material();
    const libMesh::RealTensor& rotm = mat->get_rotation_matrix();
    Tensor2Gen rotate;
    transform_tensor_format(rotm, rotate);
    strain = rotate.transpose() * strain * rotate;

    _bulk_model->set_strain(strain);

    // it wants temperature in K
    _bulk_model->set_temperature(SimulationOptions::T);

    _solve_kp();

    data.strain_seq_num = _strain_if.simulation()->get_solve_sequence_number();
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





void
KPBulkDOS::calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
    double kT, double , const Elem* , const Point& ) const
{
  double density = 0.0;
  double derivative = 0.0;
  double derivative2 = 0.0;

  const double arg_max = 150;
  const double arg_min = -50;
  const double min_dens = 1e-64;

  const double eps = 1e-7;


  _th_el_power = 0;

  for (int i = 0; i < _ref_energies.size(); ++i)
  {
    double dens, der, der2;

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
      if (result.size() > 1)
        der = dens;
      if (result.size() > 2)
        der2 = dens;
    }
    else if (arg < arg_max)
    {
      dens = TiberMath::fermidirac_half(arg);
      if (result.size() > 1)
        der = TiberMath::fermidirac_mhalf(arg);
      if (result.size() > 2)
        der2 = (TiberMath::fermidirac_mhalf(arg + eps)-TiberMath::fermidirac_mhalf(arg - eps))/(2*eps);
    }
    else
    {
      dens = 2.0 * M_2_SQRTPI / 3.0 * std::pow(arg, 1.5);
      if (result.size() > 1)
        der = M_2_SQRTPI * std::sqrt(arg);
      if (result.size() > 2)
        der2 = 0.5 * M_2_SQRTPI * std::pow(arg,-0.5);
    }

    dens *= Neff;
    der *= Neff / kT;
    der2 *= Neff / kT /kT;

    // is this needed???
    /*if (dens > min_dens)
      der *= Neff / kT;
    else
      der = 0.0;*/

    density += dens;
    derivative += der;
    derivative2 += der2;


    // calculate thermoelectric power
    double temp = kT / Constants::k_B;
    double Pth = dens / der * 1.5 / temp + Constants::k_Boltzmann * (1 - arg);

    _th_el_power += dens * Pth;

  }

  _th_el_power /= density;

  result[0] = density;
  if (result.size() > 1)
    result[1] = derivative;
  if (result.size() > 2)
    result[2] = derivative2;
}

void
KPBulkDOS::calculate_density_and_derivative(std::vector<double>& result,double Ef, double Epot,
    double kT, double ) const
{
  double dens = 0;
  result[0] = dens;
  if (result.size() > 1)
  {
    double der = 0;
    result[1] = der;
  }
  if (result.size() > 2)
  {
    double der2 = 0;
    result[2] = der2;
  }
}

