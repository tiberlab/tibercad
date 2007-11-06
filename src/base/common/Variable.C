// $Id$

#include "boost/regex.hpp"
#include "boost/algorithm/string/trim.hpp"

#include "Variable.h"
#include "Utils.h"


using namespace std;


Variable::VariableMap
Variable::_variables;


//Variable*
//Variable::get_variable(const string& var)
//{
//  Variable* sw = NULL;
//  VariableMap::iterator it(_variables.find(var));
//  if (it != _variables.end())
//    sw = (it->second).first;

//  return sw;
//}



bool
Variable::is_variable(const std::string& var)
{
  bool result = false;
  VariableMap::iterator it(_variables.find(var));
  if (it != _variables.end())
    result = true;

  return result;
}



void
Variable::set_variable_value(const string& var, double value)
{
  // we have to set the value in all objects
  VariableMap::iterator it(_variables.lower_bound(var));
  const VariableMap::iterator end(_variables.upper_bound(var));
  for (; it != end; ++it)
  {
    Variable* sw = (it->second).first;
    sw->set_variable_value(value, (it->second).second);
  }
}



double
Variable::get_variable_value(const string& var)
{
  double val = 0.0;

  // we can get the value just from the first object as it is
  // the same in all of them
  VariableMap::iterator it(_variables.find(var));
  if (it != _variables.end())
  {
    Variable* sw = (it->second).first;
    val = sw->get_variable_value((it->second).second);
  }

  return val;
}


double
Variable::check_and_register(const string& s, double defaultval, ID id)
{
  if (s.size() > 1)
  {
    if (s[0] == '@')
    {
      string str(s);
      boost::algorithm::trim(str);

      // the regexp to match
      static const boost::regex regexp("@([a-zA-Z0-9]+)(\\[([^\\]]+)\\])?");

      boost::cmatch matches;
      if (boost::regex_match(str.c_str(), matches, regexp))
      {
        int n = matches.size();
        if (n > 1)
        {
          // get the variable name and register
          string name(matches[1].first, matches[1].second);

          if (n > 3)
          {
            // get the default value
            string val(matches[3].first, matches[3].second);
            
            if (val.size() > 0)
            {
              defaultval = Utils::convert<double>(val);

              // push the new default value to all other objects using
              // the same variable
              set_variable_value(name, defaultval);
            }
            else
            {
              // we try to get the default value from an already registered
              // object, if present
              if (is_variable(name))
                defaultval = get_variable_value(name);
            }
          }


          // finally register the object
          _variables.insert(VariableMap::value_type(name,
                pair<Variable*, ID>(this, id)));
        }
      }

    }
    else
      defaultval = Utils::convert<double>(s);
  }

  return defaultval;
}
