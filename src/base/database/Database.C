// $Id$

#include "Database.h"
#include "InitFailedException.h"


#include <boost/filesystem/operations.hpp>



bool
Database::is_alloy(const std::string& name) const
{
  if (name == "AlGaAs")  return true;
  if (name == "InGaAs")  return true;
  if (name == "AlInAs")  return true;
  if (name == "AlGaN")  return true;
  if (name == "InGaN")  return true;
  if (name == "AlInN")  return true;

  return false;
}



void
Database::get_alloy_components(const std::string& alloy,
    std::string& comp_A, std::string& comp_B) const
{
  if (alloy == "AlGaAs") { comp_A = "AlAs"; comp_B = "GaAs"; }
  if (alloy == "InGaAs") { comp_A = "InAs"; comp_B = "GaAs"; }
  if (alloy == "AlInAs") { comp_A = "AlAs"; comp_B = "InAs"; }
  if (alloy == "AlGaN") { comp_A = "AlN"; comp_B = "GaN"; }
  if (alloy == "InGaN") { comp_A = "InN"; comp_B = "GaN"; }
  if (alloy == "AlInN") { comp_A = "AlN"; comp_B = "InN"; }
}


void
Database::check_data_file(const std::string& name) const
{
  using namespace boost::filesystem;

  path datafile(name);
  if (!exists(datafile) || is_directory(datafile))
  {
    std::string msg("Database: cannot find material data file ");
    msg += name;
    throw InitFailedException(msg);
  }
}
