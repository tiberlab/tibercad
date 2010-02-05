// $Id$

#include "Database.h"
#include "Utils.h"
#include "Messages.h"
#include "DatabaseException.h"

#include "getpot.h"
#include "dense_vector.h"

#include <boost/filesystem/operations.hpp>

#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;


string
Database::_default_path = "";


string
Database::_path = "";



Database::Database(const string& material,
    const string& datafile)
  : _section(""),
    _file(NULL),
    _is_alloy(false),
   _mixing_type(VCA)
{
  set_material(material, datafile);

  // check if it is an alloy and create sub-databases
  if (is_alloy())
  {
    vector<string> comp;
    get_alloy_components(comp);

    for (size_t i = 0; i < comp.size(); i++)
      _comp_db.push_back(Database(comp[i]));
  }
}


Database::Database(const Database& other)
  : _file(NULL)
{
  Database::operator=(other);
}


Database::~Database(void)
{
  close();
}


Database&
Database::operator=(const Database& rhs)
{
  if (&rhs != this)
  {
    _material = rhs._material;
    _datafile = rhs._datafile;
    _is_alloy = rhs._is_alloy;
    _mixing_type = rhs._mixing_type;
    _comp_db = rhs._comp_db;
  }

  return *this;
}


void
Database::set_material(const string& material,
    const string& datafile)
{
  string df(datafile);
  if (df.size() == 0)
    df = get_data_file(material);

  if ((_material != material) || (_datafile != df))
  {
    _material = material;
    _datafile = df;

    open();
    if (_file->have_variable("alloy"))
      _is_alloy = true;

    close();
  }
}


void
Database::do_open(void) const
{
  assert(_file == NULL);

  if (!check_data_file(_datafile))
  {
    string msg("Cannot open material data file ");
    msg += _datafile;
    throw DatabaseException(msg);
  }

  _file = new GetPot(_datafile);
}


void
Database::close()
{
  delete _file;
  _file = NULL;
}


void
Database::set_section(const string& section) const
{
  open();

  _section = section;
  if (_section.size() != 0)
    _file->set_prefix(_section + "/");
  else
    _file->set_prefix("");

  for (size_t i = 0; i < _comp_db.size(); i++)
    _comp_db[i].set_section(section);
}



void
Database::get_alloy_components(string& comp_A,
    string& comp_B) const
{
  open();
  comp_A = (*_file)("comp_A", "");
  comp_B = (*_file)("comp_B", "");
}


void
Database::get_alloy_components(vector<string>& comp) const
{
  open();
  comp.resize(2);
  comp[0] = (*_file)("comp_A", "");
  comp[1] = (*_file)("comp_B", "");
}


bool
Database::check_data_file(const string& name) const
{
  bool ans = true;

  ifstream infile;
  infile.open(name.c_str());
  if (infile.fail() || !infile.good() || (infile.rdbuf()->in_avail() == 0))
    ans = false;

  return ans;
}


void
Database::set_search_path(const string& path)
{
  if (path.size() > 0)
  {
    boost::filesystem::path p(path, boost::filesystem::native);
    if (!boost::filesystem::exists(p) || !boost::filesystem::is_directory(p))
    {
      string msg("\'");
      msg += path + "\' is not a valid directory for searchpath";
      throw DatabaseException(msg);
    }
  }

  _path = path;
}


const string&
Database::get_search_path(void)
{
  if (_path.size() != 0) return _path;

  return _default_path;
}


void
Database::set_default_search_path(const string& path)
{
  if (path.size() > 0)
  {
    boost::filesystem::path p(path, boost::filesystem::native);
    if (!boost::filesystem::exists(p) || !boost::filesystem::is_directory(p))
    {
      string msg("\'");
      msg += path + "\' is not a valid directory for default searchpath";
      throw DatabaseException(msg);
    }
  }

  _default_path = path;
}


const string
Database::get_data_file(const string& material) const
{
  string s(_path);
  s += "/" + material + ".dat";

  if ((_path.size() == 0) || !check_data_file(s))
  {
    s = _default_path + "/" + material + ".dat";

    if ((_default_path.size() == 0) || (!check_data_file(s)))
    {
      string msg("Cannot find material data file ");
      msg += material + ".dat";
      throw DatabaseException(msg);
    }
  }

  return s;
}




void
Database::require_variable(const string& variable) const
{
  if (!has_variable(variable))
  {
    string msg("Variable \'");
    msg += variable + "\' is required in section \'" + _section
      + "\' of material data file " + _datafile;
    throw DatabaseException(msg);
  }
}



template <typename T>
T
Database::get(const string& variable, T default_value,
    bool required) const
{
  open();

  T result;

  if (is_alloy() && (_mixing_type != NONE))
  {
    size_t n = get_number_of_components();
    result = _comp_db[0].get(variable, default_value, required);
    for (size_t i = 1; i < n; i++)
      if (result != _comp_db[i].get(variable, default_value, required))
      {
        ostringstream os;
        os << "In database of " << get_material() << ": parameter "
          << variable << " in section " << get_section()
          << " has different values in the alloy components.";
        Messages::warning(os.str());
      }
  }
  else
  {
    if (required) require_variable(variable);
    result = (*_file)(variable.c_str(), default_value);
  }

  return result;
}





template <typename T>
void
Database::get(const string& variable, vector<T>& data, bool required) const
{
  open();

  if (required) require_variable(variable);
  else if (!has_variable(variable)) return;

  size_t n = data.size();
  string s(get(variable, ""));
  Utils::extract_vector(s, data);

  if ((n > 0) && (data.size() != n))
  {
    ostringstream msg;
    msg << "Variable \'" << variable << "\' in section \'" << _section
      << "\' of material data file " << _datafile
      << " has to provide a vector with " << n << " components";
    throw DatabaseException(msg.str());
  }
}



template <typename T>
void
Database::get(const string& variable,
    vector<vector<T> >& data, bool required) const
{
  open();

  if (required) require_variable(variable);
  else if (!has_variable(variable)) return;

  string s(get(variable, ""));
  vector<string> vec;
  Utils::extract_vector(s, vec);

  size_t n = vec.size();
  if ((data.size() > 0) && (data.size() != n))
  { 
    ostringstream msg;
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
      ostringstream msg;
      msg << "Row " << (i + 1) << " of variable \'" << variable
        << "\' in secton \'" << _section
        << "\' of material data file " << _datafile
        << " has to have " << ns << " components";
      throw DatabaseException(msg.str());
    }
  }
}



bool
Database::has_variable(const string& variable) const
{
  open();

  return _file->have_variable(variable.c_str());
}


void
Database::set_alloy_composition(std::vector<double>& fractions)
{
  _comp_fractions = fractions;

  size_t n = fractions.size();
  if (n != get_number_of_components())
  {
    ostringstream os;
    os << "You tried to assign molar fractions for " << n << " components "
      << "to the alloy " << get_material() << " which has only "
      << get_number_of_components() << " components.";
    throw DatabaseException(os.str());
  }

  double tot = 0;
  for (size_t i = 0; i < n; i++)
  {
    tot += fractions[i];
  }
  if (!Utils::almost_equal::compare(tot, 1.0))
  {
    ostringstream os;
    os << "The molar fractions assigned "
      << "to the alloy " << get_material() << " do not sum to unity ("
      << "sum = " << tot << ").";
    throw DatabaseException(os.str());
  }
}



string
Database::get(const string& variable,
    const string& default_value, bool required) const
{
  open();
  string result;

   if (is_alloy() && (_mixing_type != NONE))
   {
     size_t n = get_number_of_components();
     result = _comp_db[0].get(variable, default_value, required);
     for (size_t i = 1; i < n; i++)
       if (result != _comp_db[i].get(variable, default_value, required))
       {
         ostringstream os;
         os << "In database of " << get_material() << ": parameter "
           << variable << " in section " << get_section()
           << " has different values in the alloy components.";
         Messages::warning(os.str());
       }
   }
   else
   {
     if (required) require_variable(variable);
     result = (*_file)(variable.c_str(), default_value);
   }

   return result;
}


string
Database::get(const string& variable,
    const char* default_value, bool required) const
{
  return get(variable, string(default_value), required);
}



//
// the specializations of get(...) for double do the mixing in case of alloys
//
// TODO for now this does only VCA
//

template <>
double
Database::get(const string& variable, double default_value,
    bool required) const
{
  open();

  double result = 0;

  if (is_alloy() && (_mixing_type != NONE))
  {
    size_t n = get_number_of_components();
    for (size_t i = 0; i < n; i++)
      result += _comp_fractions[i] * _comp_db[i].get(variable, default_value, required);

    if (n == 2)
    {
      double bow = (*_file)(string("bow_" + variable).c_str(), 0.0);
      result -= bow * _comp_fractions[0] * _comp_fractions[1];
    }

  }
  else
  {
    if (required) require_variable(variable);
    result = (*_file)(variable.c_str(), default_value);
  }

  return result;
}


template <>
void
Database::get(const string& variable, vector<double>& data, bool required) const
{
  open();

  if (is_alloy() && (_mixing_type != NONE))
  {
    size_t n = get_number_of_components();

    DenseVector<double> tmp(data);
    _comp_db[0].get(variable, tmp.get_values(), required);
    tmp.scale(_comp_fractions[0]);

    DenseVector<double> result(tmp);

    for (size_t i = 1; i < n; i++)
    {
      tmp = data;
      _comp_db[i].get(variable, tmp.get_values(), required);
      if (tmp.size() != result.size())
      {
        ostringstream os;
        os << "Array " << variable << " has different size "
          << "in the databases of the alloy components of "
          << get_material() << ".";
        throw DatabaseException(os.str());
      }
      result.add(_comp_fractions[i], tmp);
    }

    if (n == 2)
    {
      size_t nr = result.size();
      DenseVector<double> bow(nr);
      string s((*_file)(string("bow_" + variable).c_str(),""));
      Utils::extract_vector(s, bow.get_values());
      if (bow.size() == 1)
        bow.get_values() = vector<double>(nr, bow(0));
      result.add(-(_comp_fractions[0] * _comp_fractions[1]), bow);
    }

    data = result.get_values();

  }
  else
  {
    if (required) require_variable(variable);
    else if (!has_variable(variable)) return;

    size_t n = data.size();
    string s(get(variable, ""));
    Utils::extract_vector(s, data);

    if ((n > 0) && (data.size() != n))
    {
      ostringstream msg;
      msg << "Variable \'" << variable << "\' in section \'" << _section
        << "\' of material data file " << _datafile
        << " has to provide a vector with " << n << " components";
      throw DatabaseException(msg.str());
    }
  }
}



template <>
void
Database::get(const string& variable,
    vector<vector<double> >& data, bool required) const
{
  open();

  if (is_alloy() && (_mixing_type != NONE))
  {
    size_t n = get_number_of_components();

    vector<vector<double> > tmp(data);
    _comp_db[0].get(variable, tmp, required);

    vector<vector<double> > result(tmp);

    for (size_t i = 1; i < n; i++)
    {
      tmp = data;
      _comp_db[i].get(variable, tmp, required);
      if (tmp.size() != result.size())
      {
        ostringstream os;
        os << "Array " << variable << " has different size "
          << "in the databases of the alloy components of "
          << get_material() << ".";
        throw DatabaseException(os.str());
      }
      for (size_t j = 0; j < tmp.size(); ++j)
      {
        if (tmp[j].size() != result[j].size())
        {
          ostringstream os;
          os << "Array " << variable << " has different size "
            << "in the databases of the alloy components of "
            << get_material() << ".";
          throw DatabaseException(os.str());
        }
        for (size_t k =0; k < tmp[j].size(); k++)
          result[j][k] += _comp_fractions[i] * tmp[j][k];
      }
    }

    /*
    if (n == 2)
    {
      size_t nr = result.size();
      DenseVector<double> bow(nr);
      string s((*_file)(string("bow_" + variable).c_str(),""));
      Utils::extract_vector(s, bow.get_values());
      if (bow.size() == 1)
        bow.get_values() = vector<double>(nr, bow(0));
    }
    */

    data = result;
  }
  else
  {

    if (required) require_variable(variable);
    else if (!has_variable(variable)) return;

    string s(get(variable, ""));
    vector<string> vec;
    Utils::extract_vector(s, vec);

    size_t n = vec.size();
    if ((data.size() > 0) && (data.size() != n))
    {
      ostringstream msg;
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
        ostringstream msg;
        msg << "Row " << (i + 1) << " of variable \'" << variable
          << "\' in secton \'" << _section
          << "\' of material data file " << _datafile
          << " has to have " << ns << " components";
        throw DatabaseException(msg.str());
      }
    }
  }
}





// explicit instantiations

template
int Database::get(const string&, int, bool) const;

template
bool Database::get(const string&, bool, bool) const;

//template
//const char* Database::get(const string&, const char*, bool) const;


template
void Database::get(const string&, vector<int>&, bool) const;

template
void Database::get(const string&, vector<bool>&, bool) const;

template
void Database::get(const string&, vector<string>&, bool) const;



template
void Database::get(const string&,
    vector<vector<int> >&, bool) const;

template
void Database::get(const string&,
    vector<vector<bool> >&, bool) const;

template
void Database::get(const string&,
    vector<vector<string> >&, bool) const;




