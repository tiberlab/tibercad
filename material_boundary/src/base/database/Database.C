// $Id$

#include "Database.h"
#include "Utils.h"
#include "DatabaseException.h"

#include "getpot.h"

#include <boost/filesystem/operations.hpp>

#include <fstream>
#include <iostream>
#include <sstream>



std::string
Database::_default_path = "";


Database::~Database(void)
{
  delete _file;
}


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
 
    delete _file;
    _file = new GetPot(_datafile);
  }
}


void
Database::set_section(const std::string& section)
{
  _section = section;
  if (_section.size() != 0)
    _file->set_prefix(_section + "/");
  else
    _file->set_prefix("");
}



bool
Database::is_alloy(const std::string& name) const
{
  if (_file->have_variable("alloy"))
    return true;

  return false;
}



void
Database::get_alloy_components(const std::string& alloy,
    std::string& comp_A, std::string& comp_B)
{
  set_material(alloy);
  comp_A = (*_file)("comp_A", "");
  comp_B = (*_file)("comp_B", "");
}



void
Database::get_alloy_components(std::string& comp_A,
    std::string& comp_B) const
{
  comp_A = (*_file)("comp_A", "");
  comp_B = (*_file)("comp_B", "");
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


const std::string&
Database::get_search_path(void) const
{
  if (_path.size() != 0) return _path;

  return _default_path;
}


void
Database::set_default_search_path(const std::string& path)
{
  if (path.size() > 0)
  {
    boost::filesystem::path p(path, boost::filesystem::native);
    if (!boost::filesystem::exists(p) || !boost::filesystem::is_directory(p))
    {
      std::string msg("\'");
      msg += path + "\' is not a valid directory for default searchpath";
      throw DatabaseException(msg);
    }
  }

  _default_path = path;
}


const std::string
Database::get_data_file(const std::string& material) const
{
  std::string s(_path);
  s += "/" + material + ".dat";

  if ((_path.size() == 0) || !check_data_file(s))
  {
    s = _default_path + "/" + material + ".dat";

    if ((_default_path.size() == 0) || (!check_data_file(s)))
    {
      std::string msg("Cannot find material data file ");
      msg += material + ".dat";
      throw DatabaseException(msg);
    }
  }

  return s;
}




void
Database::require_variable(const std::string& variable) const
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
Database::get(const std::string& variable, std::vector<T>& data, bool required) const
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
    std::vector<std::vector<T> >& data, bool required) const
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



bool
Database::has_variable(const std::string& variable) const
{
  return _file->have_variable(variable.c_str());
}



std::string
Database::get(const std::string& variable,
    const std::string& default_value, bool required) const
{
  if (required) require_variable(variable);
  return (*_file)(variable.c_str(), default_value);
}


std::string
Database::get(const std::string& variable,
    const char* default_value, bool required) const
{
  if (required) require_variable(variable);
  return (*_file)(variable.c_str(), std::string(default_value));
}



template <typename T>
T
Database::get(const std::string& variable, T default_value,
    bool required) const
{
  if (required) require_variable(variable);
  return (*_file)(variable.c_str(), default_value);
}





// explicit instantiations

template
double Database::get(const std::string&, double, bool) const;

template
int Database::get(const std::string&, int, bool) const;

template
bool Database::get(const std::string&, bool, bool) const;

template
const char* Database::get(const std::string&, const char*, bool) const;


template
void Database::get(const std::string&, std::vector<double>&, bool) const;

template
void Database::get(const std::string&, std::vector<int>&, bool) const;

template
void Database::get(const std::string&, std::vector<bool>&, bool) const;

template
void Database::get(const std::string&, std::vector<std::string>&, bool) const;


template
void Database::get(const std::string&,
    std::vector<std::vector<double> >&, bool) const;

template
void Database::get(const std::string&,
    std::vector<std::vector<int> >&, bool) const;

template
void Database::get(const std::string&,
    std::vector<std::vector<bool> >&, bool) const;

template
void Database::get(const std::string&,
    std::vector<std::vector<std::string> >&, bool) const;




