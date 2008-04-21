// $Id$

#include "Database.h"
#include "TiberCad.h"
#include "InitFailedException.h"


#include <fstream>
#include <iostream>



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


bool
Database::check_data_file(const std::string& name) const
{
  bool ans = true;

  std::ifstream infile;
  infile.open(name.c_str());
  if (infile.fail() || !infile.good() || (infile.rdbuf()->in_avail() == 0))
    ans = false;

  return ans;
}



const std::string
Database::get_data_file(void) const
{
  std::string s(_path);
  s += "/" + _material + ".dat";

  if (!check_data_file(s))
  {
    s = TiberCad::tiberroot + "/materials/" + _material + ".dat";

    if ((TiberCad::tiberroot.size() == 0) || (!check_data_file(s)))
    {
      std::string msg("Database: cannot find material data file ");
      msg += _material + ".dat";
      throw InitFailedException(msg);
    }
  }

  return s;
}


