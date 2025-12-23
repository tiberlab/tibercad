/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file SimulationOptions.C
 * \brief tiberCAD API implementation.
 */


#include <boost/filesystem.hpp>

#include "tibercad/base/SimulationOptions.h"
#include "tibercad/base/ModelOptions.h"
#include "tibercad/base/Variable.h"
#include "tibercad/io/Messages.h"
#include "tibercad/base/InitFailedException.h"


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
  VariableValue::check_and_register(temp, temperature);

  incomplete_ionization = opts.get_option("incomplete_ionization", true);

  _verbose = opts.get_option("verbose", 1);

  scratch_path = opts.get_option("scratchpath",".");

  // create scratch directory
  path outpath(scratch_path);

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

