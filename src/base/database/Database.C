// $Id$

#include "Database.h"
#include "TiberCad.h"
#include "InitFailedException.h"

#include "getpot.h"

#include <boost/filesystem/operations.hpp>

#include <fstream>
#include <iostream>



bool
Database::is_alloy(const std::string& name) const
{
  GetPot data(get_data_file(name));
  if (data.have_variable("alloy"))
    return true;

  return false;
}



void
Database::get_alloy_components(const std::string& alloy,
    std::string& comp_A, std::string& comp_B) const
{
  GetPot data(get_data_file(alloy));
  comp_A = data("comp_A", "");
  comp_B = data("comp_B", "");
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


void
Database::set_search_path(const std::string& path)
{
  if (path.size() > 0)
  {
    boost::filesystem::path p(path);
    if (!boost::filesystem::exists(p) || !boost::filesystem::is_directory(p))
    {
      std::string msg("Database: \'");
      msg += path + "\' is not a valid directory";
      throw InitFailedException(msg);
    }
  }

  _path = path;
}




const std::string
Database::get_data_file(const std::string& material) const
{
  std::string s(_path);
  s += "/" + material + ".dat";

  if ((_path.size() == 0) || !check_data_file(s))
  {
    s = TiberCad::tiberroot + "/materials/" + material + ".dat";

    if ((TiberCad::tiberroot.size() == 0) || (!check_data_file(s)))
    {
      std::string msg("Database: cannot find material data file ");
      msg += material + ".dat";
      throw InitFailedException(msg);
    }
  }

  return s;
}


