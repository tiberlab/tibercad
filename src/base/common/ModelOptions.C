// $Id$

#include "ModelOptions.h"

#include <sstream>

template <typename T>
void
ModelOptions::set_option(const std::string& name, const T value)
{
  std::ostringstream s;
  s << value;
  
  _options[name] = s.str();
}


int
ModelOptions::get_option(const std::string& name, int default_value) const
{
  int val = default_value;

  OptionsMap::const_iterator it = _options.find(name);

  if (it != _options.end())
    val = atoi((it->second).c_str());

  return val;
}


double
ModelOptions::get_option(const std::string& name, double default_value) const
{
  double val = default_value;

  OptionsMap::const_iterator it = _options.find(name);

  if (it != _options.end())
    val = atof((it->second).c_str());

  return val;
}


bool
ModelOptions::get_option(const std::string& name, bool default_value) const
{
  bool val = default_value;

  OptionsMap::const_iterator it = _options.find(name);

  if (it != _options.end())
  {
    if (it->second == "true")
      val = true;
    else
      val = false;

    //std::istringstream s((it->second).c_str());
    //s >> val;
  }

  return val;
}


const std::string&
ModelOptions::get_option(const std::string& name,
    const std::string& default_value) const
{
  OptionsMap::const_iterator it = _options.find(name);
  
  if (it != _options.end())
    return it->second;
  else
    return default_value;
}


const char*
ModelOptions::get_option(const std::string& name,
    const char* default_value) const
{
  OptionsMap::const_iterator it = _options.find(name);
  
  if (it != _options.end())
    return it->second.c_str();
  else
    return default_value;
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

