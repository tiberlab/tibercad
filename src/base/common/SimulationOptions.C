// $Id$

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>


#include "SimulationOptions.h"
#include "ModelOptions.h"
#include "Variable.h"
#include "Messages.h"
#include "InitFailedException.h"


double
SimulationOptions::temperature = 300.0;

double&
SimulationOptions::temp = temperature;

double&
SimulationOptions::T = temperature;

bool
SimulationOptions::incomplete_ionization = false;

int
SimulationOptions::_verbose = 1;

std::string
SimulationOptions::scratch_path = ".";



int
SimulationOptions::verbose(void)
{
  return _verbose;
}


void
SimulationOptions::initialize(const ModelOptions& opts)
{
  using namespace boost::filesystem;

  std::string temp = opts.get_option("temperature", "300.0");
  Variable::check_and_register(temp, temperature);

  incomplete_ionization = opts.get_option("incomplete_ionization", true);

  _verbose = opts.get_option("verbose", 2);

  scratch_path = opts.get_option("scratchpath",".");

  // create scratch directory
  path outpath(scratch_path, native);

  if (!exists(outpath))
  {
    // we catch any error here without doing anything yet
    try {
      create_directories(outpath);
    }
    catch (...) {}
  }

  if (!(exists(outpath) && is_directory(outpath)))
  {
    std::string msg("Cannot create or use '");
    msg += outpath.string() + "' as scratch directory.";
    throw InitFailedException(msg);
  }

}

