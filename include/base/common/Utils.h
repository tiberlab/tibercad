// $Id$

#ifndef _UTILS_H_
#define _UTILS_H_

#include "InitFailedException.h"

#include "tensor.h"

#include <typeinfo>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

namespace libMesh
{
template <typename T> class VectorValue;
template <typename T> class TypeVector;
typedef VectorValue<double> RealVectorValue;
template <typename T> class TensorValue;
typedef TensorValue<double> RealTensor;
}

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
     * The commas are optional
     *
     * The vector components can contain any character sequencies which do not
     * contain the type of braces used to define the vector. When a component
     * should contain a comma, it has to be quoted with \c ", when it should
     * contain a quote, the quote has to be escaped as \verbatim \" \endverbatim
     *
     * vec.size() = 0 if input does not contain a vector.
     */
    template <typename T>
    static void extract_vector(const std::string& input, std::vector<T>& vec);


    //! Extract a real 3D vector from a string
    /*!
     * The input string has to be of the form
     * \li \verbatim ( val1, val2, val3 ) \endverbatim
     * \li \verbatim [ val1, val2, val3 ] \endverbatim
     * \li \verbatim { val1, val2, val3 } \endverbatim
     * \li \verbatim   val1, val2, val3  \endverbatim
     *
     * The commas are optional
     *
     * \note If \c input contains only one number, \c vec will be filled with
     * that number. If it contains only two numbers, the first two components
     * of \c vec will be assigned both the first of the two.
     *
     * \return \c true if successful , \c false otherwise
     */
    static void extract_vector(const std::string& input, libMesh::TypeVector<double>& vec);


    //! Extract a real tensor from a string
    /*!
     * The input string has to be of the form
     * \li \verbatim a \endverbatim diagonal tensor, isotropic case
     * \li \verbatim ( a b ) \endverbatim diagonal tensor, unisotropic case with one
     *    different principal axis
     * \li \verbatim ( a b c ) \endverbatim  diagonal tensor
     * \li \verbatim ( a b c d e f ) \endverbatim symmetric tensor with ordering
     *      \verbatim xx yy zz yz xz xy \endverbatim
     * \li \verbatim ( a b c; d e f; g h i ) \endverbatim complete tensor
     * The semicolon can be replaced by a newline.
     *
     * Any of the syntax of the  extract_vector() method can be used for
     * the row data.
     */
    static void extract_tensor(const std::string& input, libMesh::RealTensor& tensor);


    //! Scale a vector with fractional parts to all integer vector
    static void scale_to_int(Tensor1& a, double tol = 0.001);

    //! Scale 3 vectors (as tensor) with fractional parts to all integer vector
    static void scale_to_int(Tensor2Gen &a, double tol = 0.001);


    //! Tokenize a string
    /*!
     * Cuts a string into tokens using \c delimiter (default '.') and returns
     * the tokens in \c tokens
     */
    static void tokenize(const std::string& input,
        std::vector<std::string>& tokens, const char* delimiter = ".");

    //! Cuts a string into tokens using CamelRules
    static void camel_tokenize(const std::string& input,       
        std::vector<std::string>& tokens);

    //   ! Find matching strings in a vector of strings
    //static void find_matching_strings(const std::string& regex,
    //    const std::vector<std::string>& in, std::vector<std::string>& out);


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


    //! Format a time in seconds as string
    static std::string time_to_string(double seconds);


    //! Return the directory part of a filename
    static std::string dirname(const std::string& file);


    //! Return the file basename
    static std::string basename(const std::string& file);


    //! Return the file extension
    static std::string file_extension(const std::string& file);


    //! Convert a windows style path to a cygwin POSIX pathname like
    static void convert_win32_path_to_posix(std::string& path);


    //! Read a stream until the matching (closing) symbol of a grouping pair is found
    /*!
     * \param istr the stream to parse
     * \param open the opening symbol
     * \param close the closing symbol
     */
    static const std::string get_until_matching_symbol(std::istream& istr,
        char open = '(', char close = ')');


    //! Skip whitespace in a stream
    static void skip_whitespace(std::istream& istr);


    //! Trim leading and trailing set of characters from a string
    /*!
     * As default all whitespace is trimmed.
     */
    static void trim(std::string& str,
        const std::string& chars = "");


    //! A timer class, based on the times() system call
    class Timer
    {
      public:

        //! Constructor
        /*!
         * Postcondition: the timer is reset
         */
        Timer(void);
        ~Timer(void);

        //! Reset the timer
        void reset(void);

        //! Get the elapsed time as a formatted string
        std::string elapsed_string(void);

      private:

        // The start time, in seconds
        double _start;
    };


    //! Class to print loop progress in percentage
    class Progress
    {
      public:

        //! Constructor with message line and max size
        Progress(const std::string message, unsigned int max_size,
            unsigned int step_size = 1);

        ~Progress(void);
        
        //void first_message(void);

        void progress_message(unsigned int progress);
        void progress_message(void)
        { progress_message(++_progress_counter); }

      private:
      
        std::string _message;
        unsigned int _progress_size;
        unsigned int _progress_step;
        unsigned int _progress_counter;
    };


    //! Symmetric pair: <A,B>==<B,A>, is obtained by ordering a pair<A,B> such that A<=B
    template <class T > 
    class Couple : public std::pair<T, T> 
    {
      public:
     
      Couple(const T& a, const T& b): std::pair<T,T>(a,b)
      {if(a>b){(*this).first = b ; (*this).second = a; }};

      virtual ~Couple(){};

    };

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
  char* end_p;
  double value = std::strtod(val.c_str(), &end_p);
  //if (end_p != NULL)
  //  throw InitFailedException("\"" + val + "\" is not a valid double value.");

  return(value);
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
unsigned short
Utils::convert<unsigned short>(const std::string& val)
{
  return (unsigned short) convert<double>(val);
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

  if ((val == "true") || (val == "1") || (val == "TRUE")
      || (val == "yes") || (val == "y"))
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
