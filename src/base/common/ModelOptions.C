// $Id$

#include "ModelOptions.h"

#include <sstream>
#include <iostream>


using namespace std;


ModelOptions::ModelOptions(map<const string, string> options)
  : _options(options)
{
}

template <typename T>
void
ModelOptions::set_option(const string& name, const T value)
{
  ostringstream s;
  s << value;
  
  _options[name] = s.str();
}


template <typename T>
void
ModelOptions::get_option(const string& name,
    vector<T>& vec) const
{

  OptionsMap::const_iterator it = _options.find(name);

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



ModelOptions&
ModelOptions::operator+=(const ModelOptions& rhs)
{
  OptionsMap::const_iterator it = rhs._options.begin();
  const OptionsMap::const_iterator end = rhs._options.end();

  for ( ; it != end; ++it)
    _options[it->first] = it->second;

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
ModelOptions::print_all(void) const
{
  cout << "ModelOptions content: {" << endl;
  OptionsMap::const_iterator it = _options.begin();
  const OptionsMap::const_iterator end = _options.end();

  for ( ; it != end; ++it)
    cout << it->first << " -> " << it->second << endl;

  if (!_submodels.empty())
  {
    cout << "Submodel ";
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




// explicit instantiations


template
void ModelOptions::set_option<double>(const string&, const double);

template
void ModelOptions::set_option<int>(const string&, const int);

template
void ModelOptions::set_option<bool>(const string&, const bool);

template
void ModelOptions::set_option<string>(const string&,
    const string);


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

