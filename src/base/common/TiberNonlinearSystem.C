// $Id$


#include "TiberNonlinearSystem.h"

// the implemented systems/solvers
#include "TiberNonlinLS.h"
#include "TiberNonlinTR.h"
#include "TiberNonlinBR.h"
#include "TiberNonlinPetsc.h"

#include "XMonitor.h"
#include "InitFailedException.h"


#include "equation_systems.h"

#include <cassert>

using namespace std;


TiberNonlinearSystem::TiberNonlinearSystem(libMesh::EquationSystems& es,
    const string& name, const unsigned int number)
: TiberEqSystem(),
  Parent(es, name, number),
  _n_nonlin_iterations(0),
  _final_residual_norm(1e20),
  _last_step_size(1e20),
  _xmonitor(NULL)
{
  set_type(NONLINEAR);
}



TiberNonlinearSystem*
TiberNonlinearSystem::create(libMesh::EquationSystems& es,
    const std::string& sysname, const ModelOptions& options)
{
  TiberNonlinearSystem* sys = NULL;

  std::string type(options.get_name());
  if (type.empty())
    type = "linesearch";
  type = options.get_option("type", type);

  if (type == "petsc")
    sys = &(es.add_system<TiberNonlinPetsc>(sysname));
  else if (type == "linesearch")
    sys = &(es.add_system<TiberNonlinLS>(sysname));
  else if (type == "trustregion")
    sys = &(es.add_system<TiberNonlinTR>(sysname));
  else if (type == "bankrose")
    sys = &(es.add_system<TiberNonlinBR>(sysname));
  else
  {
    std::string s = "Unknown type '" +
      type + "' for nonlinear system " + sysname;
    throw InitFailedException(s);
  }

  assert(sys != NULL);
  sys->set_options(options);

  return sys;
}


void
TiberNonlinearSystem::solve(void)
{
  if (get_options().get_option("xmonitor", false))
  {
    _xmonitor = XMonitor::create(string(get_options().get_option("simulation", "?"))
        + ": nonlinear convergence monitor");
    _xmonitor->set_axis_labels("iteration nr.", "Logarithm of residual norm");
  }
  
  do_solve();
  
  
  delete _xmonitor;
  _xmonitor = NULL;
}


void
TiberNonlinearSystem::draw_point(double iteration, double error, bool logarithm)
{
  if (_xmonitor != NULL)
  {
    if (logarithm)
      error = log10(error);

    _xmonitor->draw_point(iteration, error);
  }
}
