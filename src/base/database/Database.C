// $Id$

#include "Database.h"
#include "TiberCad.h"
#include "InitFailedException.h"

#include "getpot.h"

#include <boost/filesystem/operations.hpp>

#include <fstream>
#include <iostream>




void
Database::set_material(const std::string& material,
    const std::string& datafile)
{
  std::string df(datafile);
  if (df.size() == 0)
    df = get_data_file(material);

  if ((_material != material) || (_datafile != df))
  {
    _material = material;
    _datafile = df;

    _file = GetPot(_datafile);
  }
}


void
Database::set_section(const std::string& section)
{
  _section = section;
  if (_section.size() != 0)
    _file.set_prefix(_section + "/");
  else
    _file.set_prefix("");
}



bool
Database::is_alloy(const std::string& name) const
{
  if (_file.have_variable("alloy"))
    return true;

  return false;
}



void
Database::get_alloy_components(const std::string& alloy,
    std::string& comp_A, std::string& comp_B)
{
  set_material(alloy);
  comp_A = _file("comp_A", "");
  comp_B = _file("comp_B", "");
}



void
Database::get_alloy_components(std::string& comp_A,
    std::string& comp_B) const
{
  comp_A = _file("comp_A", "");
  comp_B = _file("comp_B", "");
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
    boost::filesystem::path p(path, boost::filesystem::native);
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
      int* i = NULL;
      *i = 5;
      throw InitFailedException(msg);
    }
  }

  return s;
}




