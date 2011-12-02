// $Id$

#include "ParticleDensity.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Device.h"
#include "TiberMath.h"
#include "Embracing.h"
#include "Messages.h"
#include "Constants.h"


using namespace std;


const double
ParticleDensity::MINDENSITY = 1e-64;


ParticleDensity::ParticleDensity(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _name(""),
  _charge(-1),
  _statistics(TiberCad::BOLTZMANN),
  _use_quantum(false),
  _is_quantum(false),
  _elem(NULL),
  _density(-1.0),
  _density_derivative(-1.0),
  _gamma(1.0),
  _embracing(NULL),
  _add_continuum(true)
{

}


ParticleDensity::ParticleDensity(const string& name,
    TiberCad::Statistics statistics)
: PhysicalModelInterface(ModelOptions()),
  _name(name),
  _statistics(statistics),
  _use_quantum(false),
  _is_quantum(false),
  _elem(NULL),
  _density(-1.0),
  _density_derivative(-1.0),
  _gamma(1.0),
  _embracing(NULL),
  _add_continuum(true)
{
  if (name == "electron")
    _charge = -1;
  else if (name == "hole")
    _charge = 1;
  //else
  //  Messages::warning("The particle \'" + name + "\' is not known.");
}




ParticleDensity*
ParticleDensity::create(const ModelOptions& options)
{
  return new ParticleDensity(options);
}


void
ParticleDensity::do_init(void)
{
  _name = get_option("particle", _name);

  if (_name == "electron")
    _charge = -1;
  else if (_name == "hole")
    _charge = 1;

  _charge = get_option("charge", _charge);

  string stat("boltzmann");
  stat = get_option("statistics", stat);
  if (stat == "boltzmann")
    _statistics = TiberCad::BOLTZMANN;
  else if (stat == "fermidirac")
    _statistics = TiberCad::FERMIDIRAC;


  vector<string> qd;
  get_option("quantum_density", qd);
  for (size_t i = 0; i < qd.size(); i++)
    add_quantum_density(qd[i]);

  if (get_quantum_simulation() != NULL)
  {
    SimulationInterface* owner =
        SimulationInterface::get_simulation(get_simulator_id());

    ModelOptions::const_submodel_iterator embit(get_options().submodels_begin("embracing"));
    if (embit != get_options().submodels_end("embracing"))
    {
      Embracing* emb =
        owner->create_embracing_region(
            get_quantum_simulation(), embit->second, true);
      set_embracing(emb);
    }

    if (owner->has_environment())
    {
      Device& dev = owner->get_environment().get_device();
      dev.extract_physical_regions(get_option("barrier_regions", ""), _barrier_ids);
    }

    _add_continuum = get_option("add_continuum_in_well", _add_continuum);
  }
}


void
ParticleDensity::add_quantum_density(const std::string& name)
{
  if (name != "")
  {
    SimulationInterface* qd = SimulationInterface::find_simulation(name);
    if (qd == NULL)
    {
      string msg("ParticleDensity: ");
      msg += "no quantum density simulation '" + name + "' found.";
      throw InitFailedException(msg);
    }

    // we assume that the density variable has this name:
    string density_name("QuantumDensity");

    ID density_id = qd->get_solution_id(density_name);

    // We let it override with a more specific name
    if (_name == "electron")
      density_name = "elDensity";
    else if (_name == "hole")
      density_name = "hlDensity";

    ID spec_id = qd->get_solution_id(density_name);
    if (spec_id != INVALID_ID)
      density_id = spec_id;


    if (density_id == INVALID_ID)
    {
      string msg("ParticleDensity: ");
      msg += "quantum density simulation '" + name +
        "' has no variable '" + density_name + "'";
      throw InitFailedException(msg);
    }

    // at this point we have for sure a quantum density simulation

    // we take the eigenenergies to add a continuum
    ID cont_id = qd->get_solution_id("EigenEnergy");

    _quantum_density.push_back(qd);
    _density_ids.push_back(density_id);
    _3D_edge.push_back(cont_id);

    use_quantum_density();

  }
}




template <>
inline
void
ParticleDensity::classical_density<TiberCad::BOLTZMANN>(void)
{
  _gamma = 1.0;
  _density = _N_eff * exp(_argument);
  _density_derivative = _density / _kT;
}



template <>
inline
void
ParticleDensity::classical_density<TiberCad::FERMIDIRAC>(void)
{
  const double arg_max = 150;
  const double arg_min = -50;

  if (_argument < arg_min)
  {
    _density = exp(_argument);
    _density_derivative = _density;
  }
  else if (_argument < arg_max)
  {
    _density = TiberMath::fermidirac_half(_argument);
    _density_derivative = TiberMath::fermidirac_mhalf(_argument);
  }
  else
  {
    _density = 2.0 * M_2_SQRTPI / 3.0 * std::pow(_argument, 1.5);
    _density_derivative = M_2_SQRTPI * std::sqrt(_argument);
  }

  _density *= _N_eff;
  if (_density > 0.1 * MINDENSITY)
    _density_derivative *= _N_eff / _kT;
  else
    _density_derivative = 0.0;
}





bool
ParticleDensity::quantum_density(void)
{
  bool flag = false;
  _density = 0.0;

  double qdens = 0.0;
  // put it to something small and then check for a slightly bigger number
  double continuum = -1000.0;

  for (size_t i = 0; i < _quantum_density.size(); i++)
  {
    vector<Point> p(1, _p);
    vector<double> values(1, 0.0);

    if (_quantum_density[i]->is_solved())
      flag |= _quantum_density[i]->get_solution(_elem, _density_ids[i], values, p);

    qdens += values[0];

    if (flag && (_3D_edge[i] != INVALID_ID) && _add_continuum)
    {
      map<ID, vector<double> > tmp;
      tmp[_3D_edge[i]] = vector<double>();
      _quantum_density[i]->get_solution(tmp);
      size_t n = tmp[_3D_edge[i]].size();
      if (n > 0)
        continuum = max(continuum, tmp[_3D_edge[i]][n - 1]);
    }
  }

  if (flag)
  {
    ID subdomid = _elem->subdomain_id();

    if ((_barrier_ids.size() > 0) && _barrier_ids.count(subdomid))
    {
      // in the barrier (if specified), add the classical density
      switch (_statistics)
      {
        case TiberCad::FERMIDIRAC:
          classical_density<TiberCad::FERMIDIRAC>();
          break;
        default: // Boltzmann
          classical_density<TiberCad::BOLTZMANN>();
          break;
      }
    }
    else if (continuum > -999.0)
    {
      // in the well (or everywhere, if no barrier has been specified)
      // add a continuum from the next available energy level


      // for positive charge we have to change sign
      if (_charge > 0)
        continuum *= -1.0;

      // we do not accept it if it is lower than the bulk band edge
      continuum = max(_E, continuum);

      _argument = (_E_F - continuum) / _kT;

      switch (_statistics)
      {
        case TiberCad::FERMIDIRAC:
          classical_density<TiberCad::FERMIDIRAC>();
          break;
        default: // Boltzmann
          classical_density<TiberCad::BOLTZMANN>();
          break;
      }
      _argument = (_E_F - _E) / _kT;

    }
  }

  _density += qdens;

  return flag;
}





void
ParticleDensity::calculate_density(void)
{
  /* classical density is calculated when
   * - either _use_quantum is false
   * - or quantum_density() returns false
   */
  if (!_use_quantum || !(_is_quantum = quantum_density()))
  {
    switch (_statistics)
    {
      case TiberCad::FERMIDIRAC:
        classical_density<TiberCad::FERMIDIRAC>();
        break;
      default: // Boltzmann
        classical_density<TiberCad::BOLTZMANN>();
        break;
    }
  }
  else if (_embracing != NULL)
  {
    // we need to do a mixing
    double x = _embracing->get_mixing_coefficient(_elem, _p);
    if (x < 1.0)
    {
      double dens = x * _density;
      switch (_statistics)
      {
        case TiberCad::FERMIDIRAC:
          classical_density<TiberCad::FERMIDIRAC>();
          break;
        default: // Boltzmann
          classical_density<TiberCad::BOLTZMANN>();
          break;
      }
      dens += (1.0 - x) * _density;
      _density = dens;
    }
  }

  if (_density > 1e-9)
    _gamma = _density / (_N_eff * exp(_argument));
  else
    _gamma = 1;
}





void
ParticleDensity::set_embracing(Embracing* embracing)
{
  if (_use_quantum)
    _embracing = embracing;
}
