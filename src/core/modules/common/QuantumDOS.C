// $Id: QuantumDOS.C 3258 2012-04-03 15:10:22Z maufder $

#include "QuantumDOS.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Device.h"
#include "Messages.h"

#include "TiberModule.h"


using namespace std;



QuantumDOS::QuantumDOS(const ModelOptions& options) :
  DensityOfStates(options),
  _add_continuum(false)
{
}


void 
QuantumDOS::prepare_submodels(void)
{
  ModelOptions opts;

  ModelOptions::submodel_iterator it (get_options().submodels_begin("classical_DOS"));
  if (it != get_options().submodels_end("classical_DOS"))
  {
    opts = it->second;
  }
  else
    opts.set_option("type", "bulk_kp");

  opts.set_option("particle", get_option("particle", "electron"));

  create_submodel(_classical, "density_of_states", opts);

}

void
QuantumDOS::do_init(void)
{

  vector<string> qd;
  get_option("quantum_density", qd);
  get_option("quantum_simulation", qd);
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
    if (_add_continuum)
      throw InitFailedException("'add_continuum_in_well' is currently "
          "broken in quantum DOS");
  }

  reference_energy() = _classical->get_reference_energy();
  effective_mass() = _classical->get_effective_mass();
  effective_dos() = _classical->get_effective_dos();
}



void
QuantumDOS::do_reinit(const Elem* elem)
{
  reference_energy() = _classical->get_reference_energy();
  effective_mass() = _classical->get_effective_mass();
  effective_dos() = _classical->get_effective_dos();
}



void
QuantumDOS::add_quantum_density(const std::string& name)
{
  if (name != "")
  {
    SimulationInterface* qd = SimulationInterface::find_simulation(name);
    if (qd == NULL)
    {
      string msg("Quantum DOS: ");
      msg += "no quantum density simulation '" + name + "' found.";
      throw InitFailedException(msg);
    }

    // we assume that the density variable has this name:
    string density_name("QuantumDensity");

    ID density_id = qd->get_solution_id(density_name);

    // We let it override with a more specific name
    if (get_particle() == 'e')
      density_name = "elDensity";
    else if (get_particle() == 'h')
      density_name = "hlDensity";

    ID spec_id = qd->get_solution_id(density_name);
    if (spec_id != INVALID_ID)
      density_id = spec_id;

    if (get_particle() == 'e')
      density_name = "eDensity";
    else if (get_particle() == 'h')
      density_name = "hDensity";

    spec_id = qd->get_solution_id(density_name);
    if (spec_id != INVALID_ID)
      density_id = spec_id;



    if (density_id == INVALID_ID)
    {
      string msg("Quantum DOS: ");
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






std::pair<double, double>
QuantumDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice, const Elem* elem, const Point& p) const
{
  double density = 0.0;
  double der = 0.0;

  // true indicates that we obtained a quantum density
  bool flag = false;

  // an artificial continuum band edge (energy of next state)
  // (put it to something small and then check for a slightly bigger number)
  double continuum = -1000.0;

  if (!_quantum_density.empty() && has_quantum_density())
  {

    double qdens = 0.0;

    for (size_t i = 0; i < _quantum_density.size(); i++)
    {
      vector<Point> pt(1, p);
      vector<double> values(1, 0.0);

      if (_quantum_density[i]->is_solved())
        flag |= _quantum_density[i]->get_solution(elem, _density_ids[i], values, pt);

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

    density += qdens;
  }


  if (flag)
  {
    ID subdomid = elem->subdomain_id();

    if ((_barrier_ids.size() > 0) && _barrier_ids.count(subdomid))
    {
      pair<double, double> cl(_classical->get_occupied_density_and_derivative(
            Ef, Epot, kT, elem, p, kTlattice));
      density += cl.first;
      der += cl.second;
    }
    else if (continuum > -999.0)
    {
      // in the well (or everywhere, if no barrier has been specified)
      // add a continuum from the next available energy level

      double E = get_reference_energy()[0] + Epot;

      // we do not accept it if it is lower than the bulk band edge
      // (for positive charge we have to change sign)
      if (get_particle() == 'e')
        continuum = max(E, continuum);
      else
        continuum = min(E, continuum);


      double cl_ref = _classical->get_reference_energy()[0];
      pair<double, double> cl(_classical->get_occupied_density_and_derivative(
            Ef, continuum - cl_ref, kT, elem, p, kTlattice));
      density += cl.first;
      der += cl.second;
    }
  }
  else
  {
    //
    // here we get if there was no quantum density
    //

    pair<double, double> cl(_classical->get_occupied_density_and_derivative(
          Ef, Epot, kT, elem, p, kTlattice));
    density += cl.first;
    der += cl.second;
    
  }

  return make_pair(density, der);
}

std::pair<double, double>
QuantumDOS::calculate_density_and_derivative(double Ef, double Epot,
    double kT, double kTlattice) const
{
  return make_pair(0, 0);
}

void
QuantumDOS::set_embracing(Embracing* embracing)
{
  if (is_quantum_density())
    _embracing = embracing;
}


void
QuantumDOS::do_print_info(void)
{
  Messages::info("classical DOS: " + _classical->get_name());
  _classical->print_info();
}
