// $Id$

#include "tibercad/io/Database.h"
#include "tibercad/utils/Utils.h"
#include "tibercad/io/Messages.h"
#include "DatabaseException.h"

#include "getpot.h"
#include "dense_vector.h"

#include <boost/version.hpp>
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
    const string& datafile,
    const ModelOptions& options)
  : _section(""),
    _file(NULL),
    _options(options),
    _is_alloy(false),
    _is_interface(false),
   _mixing_type(VCA)
{
  set_material(material, datafile);

  // check if it is an alloy and create sub-databases
  if (is_alloy() || is_interface())
  {
    vector<string> comp;
    get_components(comp);

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
    _options = rhs._options;
    _is_alloy = rhs._is_alloy;
    _is_interface = rhs._is_interface;
    _mixing_type = rhs._mixing_type;
    _comp_db = rhs._comp_db;
    _comp_fractions = rhs._comp_fractions;
    //set_material(_material, _datafile); // Does not work
  }

  return *this;
}


void
Database::set_material(const string& material,
    const string& datafile)
{

  _material = material;

  if (_material.find("%", 0) != string::npos)
    _is_interface = true;
  
  if (_options.find_option("alloy"))
    _is_alloy = true;


  string df(datafile);
  if (df.size() == 0)
  {
    
    df = find_file(material + ".dat");

    // we would like to check different filenames 
    // e.g. for alloys InGaP = GaInP
    // So we try to invert the first two atom species
    if (!check_data_file(df) && !is_interface())
    {
      std::vector<std::string> elems;
      Utils::camel_tokenize(material, elems);
      if (elems.size() >= 2)
      {
        ostringstream mat;
        mat << elems[1] << elems[0];

        for (unsigned int k = 2; k < elems.size(); k++)  
          mat << elems[k];
        
        df = find_file(mat.str() + ".dat");

        // if it is found, we need to invert the alloy fractions accordingly
        if (check_data_file(df) && !is_interface())
        {
          std::swap(_comp_fractions[0], _comp_fractions[1]);
        }
      }

    }
      
  }
  else
    df = find_file(df);

  if (!check_data_file(df) && !is_interface() )
  {
    string msg("Cannot open material data file for ");
    msg += material;
    Messages::warning(msg);
    //throw DatabaseException(msg);
  }
  
  _datafile = df;
  
  // it might be an alloy
  open();
  if (_file->have_variable("alloy")) _is_alloy = true;
  
  close();
  
}




void
Database::do_open(void) const
{
  assert(_file == NULL);

  if (!_datafile.empty() && (!check_data_file(_datafile) && !is_interface()))
  {
    string msg("Cannot open material data file '");
    msg += _datafile + "'";
    Messages::warning(msg);
    //throw DatabaseException(msg);
  }

  if (_datafile.empty())
    _file = new GetPot();
  else
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
Database::get_components(vector<string>& comp) const
{
  if (is_alloy())
  {
    open();
    // TODO extend to any type of alloy
    comp.resize(2);
    comp[0] = (*_file)("comp_A", "");
    comp[1] = (*_file)("comp_B", "");
  }
  else if (is_interface())
  {
    comp.resize(2);
    size_t pos = _material.find("%", 0);
    comp[0] = _material.substr(0, pos);
    comp[1] = _material.substr(pos + 1, _material.size() - pos);
  }
}


bool
Database::check_data_file(const string& name)
{
  bool ans = true;

  ifstream infile;
  infile.open(name.c_str());
  //if (infile.fail() || !infile.good() || (infile.rdbuf()->in_avail() == 0))
  if (infile.fail() || !infile.good())
    ans = false;

  return ans;
}


void
Database::set_search_path(const string& path)
{
  if (path.size() > 0)
  {
#if BOOST_VERSION >= 104700
    boost::filesystem::path p(path);
#else
    boost::filesystem::path p(path, boost::filesystem::native);
#endif
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
#if BOOST_VERSION >= 104700
    boost::filesystem::path p(path);
#else
    boost::filesystem::path p(path, boost::filesystem::native);
#endif
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
Database::find_file(const string& file) const
{

  

  // Needs boost >1.41
  boost::filesystem::path p(file);
  if (p.is_relative())
    p = _path / p;

  if ((_path.size() == 0) || !check_data_file(p.string()))
  {
    p = boost::filesystem::path(_default_path) / file;

    if ((_default_path.size() == 0) || (!check_data_file(p.string())))
    {
      p = "";
    }
  }

  return p.string();
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
  data.resize(n);
  for (size_t i = 0; i < n; i++)
  {
    size_t ns = data[i].size();
    Utils::extract_vector(vec[i], data[i]);
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
Database::get(const string& variable, const string& default_value, bool required) const
{
  open();
  string result;

   if (is_alloy() && (_mixing_type != NONE))
   {
     //for alloys we always get the string of the first component
     result = _comp_db[0].get(variable, default_value, required);

    // override from alloy DB
    result = (*_file)(variable.c_str(), result);
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
    {
      result += _comp_fractions[i] * _comp_db[i].get(variable, default_value, required);
    }

    if (n == 2)
    {
      double bow = (*_file)(string("bow_" + variable).c_str(), 0.0);
      result -= bow * _comp_fractions[0] * _comp_fractions[1];
    }

    // override from alloy DB
    result = (*_file)(variable.c_str(), result);

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

    libMesh::DenseVector<double> tmp(data);
    _comp_db[0].get(variable, tmp.get_values(), required);
    tmp.scale(_comp_fractions[0]);

    libMesh::DenseVector<double> result(tmp);

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
      libMesh::DenseVector<double> bow(nr);
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
  }

  // override from alloy file if present
  if (has_variable(variable))
  {

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
  }

  if (has_variable(variable))
  {

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



void
Database::get(const std::string& variable, libMesh::RealVectorValue& data,
    bool required) const
{

  vector<double> data_vec;

  get(variable, data_vec, required);

  switch (data_vec.size())
  {
    case 1:
      data(0) = data(1) = data(2) = data_vec[0];
      break;

    case 2:
      data(0) = data(1) = data_vec[0];
      data(2) = data_vec[1];
      break;

    case 3:
      data(0) = data_vec[0];
      data(1) = data_vec[1];
      data(2) = data_vec[2];
      break;

    default:
      break;
  }

}

void
Database::get(const string& variable, libMesh::RealTensor& tensor, 
    bool required) const
{
  open();

  if (is_alloy() && (_mixing_type != NONE))
  {
    size_t n = get_number_of_components();

    libMesh::RealTensor tmp;

    for (size_t i = 0; i < n; i++)
    {
      _comp_db[i].get(variable, tmp, required);

      tensor += _comp_fractions[i] * tmp;
    }

    if (n == 2)
    {
      libMesh::RealTensor bow;
      string s((*_file)(string("bow_" + variable).c_str(),"0.0"));

      Utils::extract_tensor(s, bow);

      tensor -= bow * _comp_fractions[0] * _comp_fractions[1];

    }

  }
  else
  {
    if (required) require_variable(variable);
  }
  if (has_variable(variable))
  {
    string s(get(variable, ""));
    Utils::extract_tensor(s, tensor);
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




