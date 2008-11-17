// $Id$

#include "Database.h"
#include "Utils.h"
#include "TiberCad.h"
#include "DatabaseException.h"

#include "getpot.h"

#include <boost/filesystem/operations.hpp>

#include <fstream>
#include <iostream>
#include <sstream>




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
      std::string msg("\'");
      msg += path + "\' is not a valid directory for searchpath";
      throw DatabaseException(msg);
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
      std::string msg("Cannot find material data file ");
      msg += material + ".dat";
      throw DatabaseException(msg);
    }
  }

  return s;
}




void
Database::require_variable(const std::string& variable)
{
  if (!has_variable(variable))
  {
    std::string msg("Variable \'");
    msg += variable + "\' is required in section \'" + _section
      + "\' of material data file " + _datafile;
    throw DatabaseException(msg);
  }
}




template <typename T>
void
Database::get(const std::string& variable, std::vector<T>& data, bool required)
{
  if (required) require_variable(variable);
  else if (!has_variable(variable)) return;

  int n = data.size();
  std::string s(get(variable, ""));
  Utils::extract_vector(s, data);

  if ((n > 0) && (data.size() != n))
  {
    std::ostringstream msg;
    msg << "Variable \'" << variable << "\' in section \'" << _section
      << "\' of material data file " << _datafile
      << " has to provide a vector with " << n << " components";
    throw DatabaseException(msg.str());
  }
}



template <typename T>
void
Database::get(const std::string& variable,
    std::vector<std::vector<T> >& data, bool required)
{
  if (required) require_variable(variable);
  else if (!has_variable(variable)) return;

  std::string s(get(variable, ""));
  std::vector<std::string> vec;
  Utils::extract_vector(s, vec);

  size_t n = vec.size();
  if ((data.size() > 0) && (data.size() != n))
  { 
    std::ostringstream msg;
    msg << "Variable \'" << variable << "\' in section \'" << _section
      << "\' of material data file " << _datafile
      << " has to provide an array with " << data.size() << " rows";
    throw DatabaseException(msg.str());
  }

  data.resize(n);
  for (size_t i = 0; i < n; i++)
  {
    size_t ns = data[i].size();
    Utils::extract_vector(vec[i], data[i]);
    if ((ns > 0) && (data[i].size() != ns))
    {
      std::ostringstream msg;
      msg << "Row " << (i + 1) << " of variable \'" << variable
        << "\' in secton \'" << _section
        << "\' of material data file " << _datafile
        << " has to have " << ns << " components";
      throw DatabaseException(msg.str());
    }
  }
}





// explicit instantiations
template
void Database::get(const std::string&, std::vector<double>&, bool);

template
void Database::get(const std::string&, std::vector<int>&, bool);

template
void Database::get(const std::string&, std::vector<bool>&, bool);

template
void Database::get(const std::string&, std::vector<std::string>&, bool);


template
void Database::get(const std::string&,
    std::vector<std::vector<double> >&, bool);

template
void Database::get(const std::string&,
    std::vector<std::vector<int> >&, bool);

template
void Database::get(const std::string&,
    std::vector<std::vector<bool> >&, bool);

template
void Database::get(const std::string&,
    std::vector<std::vector<std::string> >&, bool);




