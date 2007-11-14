// $Id$

#ifndef _UTILS_H_
#define _UTILS_H_

#include <typeinfo>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

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
     * \li \verbatim   val1, val2, ..., valN  \endverbatim
     *
     * The vector components can contain any character sequencies which do not
     * contain the type of braces used to define the vector. When a component
     * should contain a comma, it has to be quoted with \c ", when it should
     * contain a quote, the quote has to be escaped as \verbatim \" \endverbatim
     */
    template <typename T>
    static void extract_vector(const std::string& input, std::vector<T>& vec);


    //! Tokenize a string
    /*!
     * Cuts a string into tokens using \c delimiter (default '.') and returns
     * the tokens in \c tokens
     */
    static void tokenize(const std::string& input,
        std::vector<std::string>& tokens, const char* delimiter = ".");


    //! A functor that checks if two double values are almost equal
    /*!
     * \em almost \em equal means the following:
     * \f[ \vert a - b \vert < \epsilon(1 + \vert a\vert) \f]
     */
    class almost_equal
    {

      public:

        //! Constructor
        /*!
         * \param eps the precision for the comparison
         *
         * Default for the precision is 1e-12
         */
        explicit almost_equal(double eps = 1e-12) : _eps(eps) {};


        typedef double first_argument_type;
        typedef double second_argument_type;
        typedef bool result_type;
        
        //! compare \c a with \c b
        bool operator()(double a, double b) const;

        //! compare \c a with \c b using precision \c eps
        bool operator()(double a, double b, double eps) const;

        //! A static version
        static bool compare(double a, double b, double eps = 1e-12);

      private:

        //! The precision
        double _eps;
    };



    //! Calculate the Bernoulli function of x
    /*!
     * Calculates \f$\frac{x}{e^x - 1}\f$
     */
    static double bernoulli(double x);



    //! Calculate the inverse of the Bernoulli function of x
    /*!
     * Calculates \f$\frac{e^x - 1}{x}\f$
     */
    static double bernoulli_inv(double x);






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
unsigned int
Utils::convert<unsigned int>(const std::string& val)
{
  return (unsigned int) convert<double>(val);
}


template<>
inline
short
Utils::convert<short>(const std::string& val)
{
  return (short) convert<double>(val);
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

  if ((val == "true") || (val == "1") || (val == "TRUE"))
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


inline
bool
Utils::almost_equal::operator()(double a, double b) const
{
  double diff = std::fabs(a - b);
  bool result = std::isless(diff, _eps * (1.0 + std::fabs(a)));
  return result;
}
  

inline
bool
Utils::almost_equal::operator()(double a, double b, double eps) const
{
  double diff = std::fabs(a - b);
  bool result = std::isless(diff, eps * (1.0 + std::fabs(a)));
  return result;
}
  

inline
bool
Utils::almost_equal::compare(double a, double b, double eps)
{
  double diff = std::fabs(a - b);
  bool result = std::isless(diff, eps * (1.0 + std::fabs(a)));
  return result;
}
  


inline
double
Utils::bernoulli(double x)
{
  double res = 1.0;

  if (std::abs(x) > 1e-12)
    res = x / (std::exp(x) - 1);

  return res;
}




inline
double
Utils::bernoulli_inv(double x)
{
  double res = 1.0;

  if (std::abs(x) > 1e-12)
    res = (std::exp(x) - 1) / x;

  return res;
}


#endif // _UTILS_H_
