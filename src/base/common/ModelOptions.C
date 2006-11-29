// $Id$

#include "ModelOptions.h"

#include <sstream>
#include <iostream>


template <typename T>
void
ModelOptions::set_option(const std::string& name, const T value)
{
  std::ostringstream s;
  s << value;
  
  _options[name] = s.str();
}


template <typename T>
void
ModelOptions::get_option(const std::string& name,
    std::vector<T>& vec) const
{

  OptionsMap::const_iterator it = _options.find(name);

  if (it != _options.end())
    Utils::extract_vector(it->second, vec);
}


template <typename T>
void
ModelOptions::get_option(const std::string& name,
    std::vector<std::vector<T> >& array) const
{
  std::vector<std::string> vec;
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





// explicit instantiations


template
void ModelOptions::set_option<double>(const std::string&, const double);

template
void ModelOptions::set_option<int>(const std::string&, const int);

template
void ModelOptions::set_option<bool>(const std::string&, const bool);

template
void ModelOptions::set_option<std::string>(const std::string&,
    const std::string);


template
void
ModelOptions::get_option<double>(const std::string& name,
    std::vector<double>& vec) const;

template
void
ModelOptions::get_option<int>(const std::string& name,
    std::vector<int>& vec) const;

template
void
ModelOptions::get_option<bool>(const std::string& name,
    std::vector<bool>& vec) const;

template
void
ModelOptions::get_option<char>(const std::string& name,
    std::vector<char>& vec) const;

template
void
ModelOptions::get_option<std::string>(const std::string& name,
    std::vector<std::string>& vec) const;


template
void
ModelOptions::get_option<double>(const std::string& name,
    std::vector<std::vector<double> >& vec) const;

template
void
ModelOptions::get_option<int>(const std::string& name,
    std::vector<std::vector<int> >& vec) const;

template
void
ModelOptions::get_option<bool>(const std::string& name,
    std::vector<std::vector<bool> >& vec) const;

template
void
ModelOptions::get_option<char>(const std::string& name,
    std::vector<std::vector<char> >& vec) const;

template
void
ModelOptions::get_option<std::string>(const std::string& name,
    std::vector<std::vector<std::string> >& vec) const;

