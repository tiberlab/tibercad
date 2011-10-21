#include "Vff.h"
#include "Messages.h"
#include "AtomisticStructure.h"

TIBER_MODULE(Vff, MODULE_NAME)


Vff*
Vff::_this = NULL;


Vff::Options::Options(void)
: boundary_conditions("free_standing"),
  substrate_plane("z"),
  boundary_tol(1.0)
{
}

Vff::Vff(const ModelOptions& options) :
  SimulationInterface(options)
{
  // there's nothing to be done
}


Vff::~Vff(void)
{
  // there's nothing to be done
}


Vff*
Vff::create(const ModelOptions& options)
{
  return new Vff(options);
}

void
Vff::parse_options()
{
  const ModelOptions& opts = get_options();
  Options& myopts = get_my_options();

  myopts.boundary_conditions = opts.get_option("boundary_conditions", "substrate");
  myopts.substrate_plane = opts.get_option("substrate_plane", "z");
  myopts.boundary_tol = opts.get_option("substrate_plane", 1.0);

}


void
Vff::do_init()
{
  Messages::info("Initializing VFF module");

  if (get_atomistic_structure()==NULL)
    throw InitFailedException("VFF: could not find atomistic structure");

  Messages::debug("Setting boundary conditions");

  set_boundary();

}


void
Vff::set_boundary(void)
{
  int n_atoms(get_atomistic_structure()->get_N_atoms());

  std::vector<unsigned int>& free_atoms = get_free_atoms();
  std::vector<double>& dof = get_dof();

  free_atoms.reserve(n_atoms);
  dof.reserve(n_atoms * 3);

  //Assign free atoms indexes according to selected boundary conditions
  if (get_my_options().boundary_conditions == "free_standing")
    {

      for (unsigned int i = 0; i < n_atoms; i++) free_atoms.push_back(i);

    }

}
