// $Id$

#include "ModelOptions.h"
#include "Messages.h"
#include "InitFailedException.h"

#include <sstream>
#include <iostream>
#include <list>
#include <algorithm>
#include <cassert>


using namespace std;


ModelOptions::ModelOptions(map<const string, string> options)
  : _options(options)
{
}


ModelOptions::ModelOptions(const ModelOptions& other)
{
  operator+=(other);

  _key = other._key;
  _name = other._name;
}



void
ModelOptions::set_name(const std::string& name)
{
  _name = name;
  set_option("name", get_option("name", name));
}


template <typename T>
void
ModelOptions::set_option(const string& name, const T& value)
{
  ostringstream s;
  s << value;

  _options[name] = s.str();
}



template <typename T>
void
ModelOptions::set_option(const string& name, const vector<T>& value)
{
  ostringstream s;
  unsigned int n = value.size();
  if (n > 0)
  {
    s << "(" << value[0];
    for (unsigned int i = 1; i < n; i++)
      s << "," << value[i];
    s << ")";
  }

  _options[name] = s.str();
}



ModelOptions::OptionsMap::const_iterator
ModelOptions::_find(const std::string& name) const
{
  OptionsMap::const_iterator it(_options.find(name));
  _used.insert(name);

  return it;
/*
  size_t idx = 0;
  for ( ; (idx < name.size()) && (name[idx] != '!'); idx++);

  // this is the most usual case, I guess
  if (idx == name.size())
    return _options.find(name);

  string s(name.substr(0, idx));

  const OptionsMap::const_iterator end(_options.end());
  OptionsMap::const_iterator it(end);

  idx++;
  for ( ; (it == end) && (idx < name.size()); idx++)
  {
    s.push_back(name[idx]);
    it = _options.find(s);
  }
  return it;
*/
}


void
ModelOptions::get_option(const string& name,
    RealVectorValue& vec) const
{
  OptionsMap::const_iterator it(_find(name));

  if (it != _options.end())
    Utils::extract_vector(it->second, vec);
}


template <typename T>
void
ModelOptions::get_option(const string& name,
    vector<T>& vec) const
{
  OptionsMap::const_iterator it(_find(name));

  if (it != _options.end())
    Utils::extract_vector(it->second, vec);
}


template <typename T>
void
ModelOptions::get_option(const string& name,
    vector<vector<T> >& array) const
{
  vector<string> vec;
  get_option(name, vec);

  int n = vec.size();
  array.resize(n);

  for (int i = 0; i < n; i++)
    Utils::extract_vector(vec[i], array[i]);

}




void
ModelOptions::add_submodel(const string& name,
    const ModelOptions& options)
{
  SubmodelMap::iterator it(_submodels.insert(SubmodelMap::value_type(name, options)));
  assert(it != _submodels.end());
  it->second.set_key(name);
}




ModelOptions&
ModelOptions::operator+=(const ModelOptions& rhs)
{
  OptionsMap::const_iterator it(rhs._options.begin());
  const OptionsMap::const_iterator end(rhs._options.end());

  for ( ; it != end; ++it)
    _options[it->first] = it->second;


  const_submodel_iterator subit(rhs.submodels_begin());
  const const_submodel_iterator subend(rhs.submodels_end());

  for ( ; subit != subend; ++subit)
    add_submodel(subit->first, subit->second);

  return *this;
}


ModelOptions&
ModelOptions::operator+=(const map<const string, string>& rhs)
{
  OptionsMap::const_iterator it = rhs.begin();
  const OptionsMap::const_iterator end = rhs.end();

  for ( ; it != end; ++it)
    _options[it->first] = it->second;


  return *this;
}



void
ModelOptions::clear(void)
{
  _options.clear();
  _submodels.clear();
}



void
ModelOptions::delete_all_submodels(void)
{
  _submodels.clear();
}



void
ModelOptions::print_all(void) const
{
  cout << "ModelOptions content: {" << endl;
  OptionsMap::const_iterator it = _options.begin();
  const OptionsMap::const_iterator end = _options.end();

  for ( ; it != end; ++it)
    cout << it->first << " -> " << it->second << endl;

  if (!_submodels.empty())
  {
    cout << "Submodels:" << endl;
    const_submodel_iterator it(_submodels.begin());
    const const_submodel_iterator end(_submodels.end());
    for ( ; it != end; ++it)
    {
      cout << it->first << ":\n";
      (it->second).print_all();
    }
  }
  cout << "}" << endl;
}



void
ModelOptions::delete_submodel(submodel_iterator it)
{
  _submodels.erase(it);
}




void
ModelOptions::delete_submodels(const std::string& name)
{
  _submodels.erase(name);
}




bool
ModelOptions::check_unused(int mode) const
{
  bool found = false;

  // shortcut to disable all tests:
  //mode = 0;

  list<std::string> unused;

  OptionsMap::const_iterator it(_options.begin());
  OptionsMap::const_iterator end(_options.end());
  for ( ; it != end; ++it)
    if (_used.find(it->first) == _used.end())
      unused.push_back(it->first);

  if (!unused.empty())
  {
    found = true;

    if (mode > 0)
    {
      ostringstream os;
      os << "found unused options in block "
          << get_key();
      if (!get_name().empty())
        os << " (" << get_name() << ")";
      os << ": ";
      list<std::string>::iterator uit(unused.begin());
      os << *uit;
      ++uit;
      const list<std::string>::iterator uend(unused.end());
      for ( ; uit != uend; ++uit)
        os << ", " << *uit;

      switch (mode)
      {
        case 1:
          Messages::warning(os.str());
          break;

        default:
          throw InitFailedException(os.str());
          break;
      }

    }
  }

  return found;
}


// explicit instantiations


template
void ModelOptions::set_option<double>(const string&, const double&);

template
void ModelOptions::set_option<int>(const string&, const int&);

template
void ModelOptions::set_option<unsigned int>(const string&, const unsigned int&);

template
void ModelOptions::set_option<bool>(const string&, const bool&);

template
void ModelOptions::set_option<string>(const string&, const string&);


template
void ModelOptions::set_option<double>(const string&, const vector<double>&);

template
void ModelOptions::set_option<int>(const string&, const vector<int>&);

template
void ModelOptions::set_option<unsigned int>(const string&,
    const vector<unsigned int>&);

template
void ModelOptions::set_option<bool>(const string&, const vector<bool>&);

template
void ModelOptions::set_option<string>(const string&, const vector<string>&);



template
void
ModelOptions::get_option<double>(const string& name,
    vector<double>& vec) const;

template
void
ModelOptions::get_option<int>(const string& name,
    vector<int>& vec) const;

template
void
ModelOptions::get_option<unsigned int>(const string& name,
    vector<unsigned int>& vec) const;

template
void
ModelOptions::get_option<short>(const string& name,
    vector<short>& vec) const;


template
void
ModelOptions::get_option<bool>(const string& name,
    vector<bool>& vec) const;

template
void
ModelOptions::get_option<char>(const string& name,
    vector<char>& vec) const;

template
void
ModelOptions::get_option<string>(const string& name,
    vector<string>& vec) const;


template
void
ModelOptions::get_option<double>(const string& name,
    vector<vector<double> >& vec) const;

template
void
ModelOptions::get_option<int>(const string& name,
    vector<vector<int> >& vec) const;

template
void
ModelOptions::get_option<bool>(const string& name,
    vector<vector<bool> >& vec) const;

template
void
ModelOptions::get_option<char>(const string& name,
    vector<vector<char> >& vec) const;

template
void
ModelOptions::get_option<string>(const string& name,
    vector<vector<string> >& vec) const;

