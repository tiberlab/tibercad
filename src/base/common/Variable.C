// $Id$

#include "tibercad/base/Variable.h"
#include "tibercad/base/TypedVariable.h"
#include "tibercad/utils/Utils.h"
#include "tibercad/base/InitFailedException.h"
#include "tibercad/base/RuntimeException.h"


#include "boost/regex.hpp"
#include "boost/algorithm/string/trim.hpp"

using namespace std;


VariableValue::VariableValueMap
VariableValue::_variables;


VariableValue::VariableValue(const std::string& name)
  : _name(name)
{
  // register the object
  _variables[name] = this;
}



VariableValue::~VariableValue(void)
{

}



void
VariableValue::clear_all(void)
{
  VariableValueMap::iterator it(_variables.begin());
  for ( ; it != _variables.end(); ++it)
    delete it->second;

  _variables.clear();
}


void
VariableValue::check_name(const std::string& name)
{
  if ((name.size() < 2) || (name[0] != '$'))
  {
    throw RuntimeException("\'" + name + "\' is not a valid variable name. "
        "VariableValues have to start with \'$\'");
  }
}


bool
VariableValue::is_variable(const std::string& var)
{
  bool result = false;
  check_name(var);
  VariableValueMap::iterator it(_variables.find(var));
  if (it != _variables.end())
    result = true;

  return result;
}



template <typename T>
void
VariableValue::set_variable_value(const string& var, const T& value)
{
  check_name(var);

  // we have to set the value in all objects
  VariableValueMap::iterator it(_variables.find(var));
  if (it != _variables.end())
  {
    TypedVariable<T>* sw =
      dynamic_cast<TypedVariable<T>*>(it->second);
    if (sw == NULL)
      throw RuntimeException("Type mismatch for variable\'"
          + var + "\'");
    sw->set_value(value);
  }
}



template <typename T>
T
VariableValue::get_variable_value(const string& var)
{
  check_name(var);

  T val;

  VariableValueMap::iterator it(_variables.find(var));
  if (it != _variables.end())
  {
    const TypedVariable<T>* sw =
      dynamic_cast<const TypedVariable<T>*>(it->second);
    if (sw == NULL)
      throw RuntimeException("Type mismatch for variable\'"
          + var + "\'");
    val = sw->get_value();
  }

  return val;
}





template <typename T>
string
VariableValue::check_and_register(const string& s, T& variable,
    const TiberModelObject* ct, InitializerBase<T>* initfunc)
{
  string name("");
  if (s.size() >= 1)
  {
    if (s[0] == '$')
    {
      string str(s);
      boost::algorithm::trim(str);

      // the regexp to match
      static const boost::regex regexp("(\\$[a-zA-Z0-9_]+)(\\[([^\\]]+)\\])?");

      boost::cmatch matches;
      if (boost::regex_match(str.c_str(), matches, regexp))
      {
        int n = matches.size();
        if (n > 1)
        {
          // get the variable name and register
          name = string(matches[1].first, matches[1].second);

          if (n > 3)
          {
            // get the default value
            string val(matches[3].first, matches[3].second);

            if (val.size() > 0)
            {
              variable = Utils::convert<T>(val);

              // push the new default value to all other objects using
              // the same variable
              set_variable_value(name, variable);
            }
            else
            {
              // we try to get the default value from an already registered
              // object, if present
              if (is_variable(name))
                variable = get_variable_value<T>(name);
            }
          }

          // does it already exist?
          TypedVariable<T>* var;
          VariableValueMap::iterator it(_variables.find(name));
          if (it != _variables.end())
          {
            VariableValue* v = it->second;
            var = dynamic_cast<TypedVariable<T>*>(v);
          }
          else
            var = new TypedVariable<T>(name);

          if (var == NULL)
            throw InitFailedException("Could not create variable \'"
                + name +"\'");

          var->register_variable(variable, ct, initfunc);

        }
      }
      else
        throw InitFailedException(s + " is not a valid variable declaration.");

    }
    else
    {
      variable = Utils::convert<T>(s);

      if (initfunc != NULL)
        (*initfunc)(variable);
    }
  }

  return(name);
}


void
VariableValue::unregister(const TiberModelObject* ct)
{
  VariableValueMap::iterator it(_variables.begin());
  for ( ; it != _variables.end(); ++it)
    (it->second)->do_unregister(ct);
}



//
// explicit instantiations
//

template
string
VariableValue::check_and_register<double>(const string&, double&,
    const TiberModelObject*, InitializerBase<double>* initfunc);


template
string
VariableValue::check_and_register<int>(const string&, int&,
    const TiberModelObject*, InitializerBase<int>* initfunc);


template
string
VariableValue::check_and_register<string>(const string&, string&,
    const TiberModelObject*, InitializerBase<string>* initfunc);


template
string
VariableValue::check_and_register<unsigned int>(const string&, unsigned int&,
    const TiberModelObject*, InitializerBase<unsigned int>* initfunc);


template
string
VariableValue::check_and_register<char>(const string&, char&,
    const TiberModelObject*, InitializerBase<char>* initfunc);


template
string
VariableValue::check_and_register<bool>(const string&, bool&,
    const TiberModelObject*, InitializerBase<bool>* initfunc);




template
void
VariableValue::set_variable_value<double>(const string&, const double&);

template
void
VariableValue::set_variable_value<int>(const string&, const int&);

template
void
VariableValue::set_variable_value<string>(const string&, const string&);

template
void
VariableValue::set_variable_value<char>(const string&, const char&);

template
void
VariableValue::set_variable_value<unsigned int>(const string&, const unsigned int&);

template
void
VariableValue::set_variable_value<bool>(const string&, const bool&);




template
double
VariableValue::get_variable_value<double>(const string&);

template
int
VariableValue::get_variable_value<int>(const string&);

template
string
VariableValue::get_variable_value<string>(const string&);

template
char
VariableValue::get_variable_value<char>(const string&);

template
unsigned int
VariableValue::get_variable_value<unsigned int>(const string&);

template
bool
VariableValue::get_variable_value<bool>(const string&);
