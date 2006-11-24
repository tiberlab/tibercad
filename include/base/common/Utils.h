// $Id$

#ifndef _UTILS_H_
#define _UTILS_H_

#include <typeinfo>
#include <string>

//! This class contains useful methods for different tasks
class Utils
{

  public:

    //! Extracts the human readable class name from a type_info object
    static std::string extract_typename(const std::type_info& info);

  private:

    //! Not to be instantiated
    Utils(void);

};


inline
std::string
Utils::extract_typename(const std::type_info& info)
{
  const char* s = info.name();

  if (s[0] == 'P')
    s += 1;

  return s;
}


#endif // _UTILS_H_
