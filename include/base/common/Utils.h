// $Id$

#ifndef _UTILS_H_
#define _UTILS_H_

#include <typeinfo>
#include <string>
#include <vector>

//! This class contains useful methods for different tasks
class Utils
{

  public:

    //! Extract the human readable class name from a type_info object
    static std::string extract_typename(const std::type_info& info);

    //! Convert from a string to another type
    template <typename T>
    static T convert(const std::string& val);

    //! Extract a vector of type \c T from a string
    /*!
     * \param input the string to parse
     * \param vec the extracted vector
     * 
     * The input string has to be of one of the following forms:
     * \li \verbatim ( val1, val2, ..., valN ) \endverbatim
     * \li \verbatim [ val1, val2, ..., valN ] \endverbatim
     * \li \verbatim { val1, val2, ..., valN } \endverbatim
     *
     * The vector components can contain any character sequencies which do not
     * contain the type of braces used to define the vector. When a component
     * should contain a comma, it has to be quoted with \c ", when it should
     * contain a quote, the quote has to be escaped as \verbatim \" \endverbatim
     */
    template <typename T>
    static void extract_vector(const std::string& input, std::vector<T>& vec);


  private:

    //! Not to be instantiated
    Utils(void);

};


//
// inline methods
//

template<>
inline
double
Utils::convert<double>(const std::string& val)
{
  return atof(val.c_str());
}


template<>
inline
int
Utils::convert<int>(const std::string& val)
{
  return (int) convert<double>(val);
}

template<>
inline
char
Utils::convert<char>(const std::string& val)
{
  return val.c_str()[0];
}


template<>
inline
bool
Utils::convert<bool>(const std::string& val)
{
  bool res;

  if ((val == "true") || (val == "1"))
    res = true;
  else
    res = false;

  return res;
}


template<>
inline
const std::string&
Utils::convert<const std::string&>(const std::string& val)
{
  return val;
}


template<>
inline
std::string
Utils::convert<std::string>(const std::string& val)
{
  return val;
}


template<>
inline
const char*
Utils::convert<const char*>(const std::string& val)
{
  return val.c_str();
}




#endif // _UTILS_H_
