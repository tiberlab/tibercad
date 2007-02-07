// $Id$

#ifndef _MODELOPTIONS_H_
#define _MODELOPTIONS_H_

#include "Utils.h"

#include <map>
#include <vector>
#include <string>

//! A class to store model options
/*!
 * All options are stored internally as strings. They are accessed by
 * providing a key (which is also a string).
 * The key is assumed to be unique.
 */
class ModelOptions
{

  public:

    //! The default constructor
    ModelOptions(void) {};

    //! Constructor which takes a map as argument
    ModelOptions(std::map<const std::string, std::string> options);

    //! The destructor
    ~ModelOptions(void) {};

    //! Check if it is empty
    bool is_empty(void);

    //! Get the value of an option
    /*!
     * \param name the name of the option
     * \param default_value the default value, which also defines
     * the type of the option
     * \return the value
     */
    template <typename T>
    T get_option(const std::string& name, T default_value) const;

    //! Get an option which is a vector of values (of the same type)
    /*!
     * \param name the name of the option
     * \param vec the vector, where the values will be stored. \c vec can
     * contain default values, but it's size will be changed according to
     * the vector found in the options.
     */
    template <typename T>
    void get_option(const std::string& name, std::vector<T>& vec) const;

    //! Get an option which is a vector of vectors (with the same type)
    /*!
     * \param name the name of the option
     * \param array the arry, where the values will be stored. \c array can
     * contain default values, but it's size will be changed according to
     * the array found in the options.
     *
     * The subvectors have to be quoted and have to use another type of
     * braces. E.g.:
     * \code {(1), "(1,2,3)", "(3,4)"} \endcode
     */
    template <typename T>
    void get_option(const std::string& name,
        std::vector<std::vector<T> >& vec) const;

    //! Check if an option is present
    bool find_option(const std::string& name) const;

    //! Set an option
    /*!
     * This method is probably of no use...
     */
    template <typename T>
    void set_option(const std::string& name, const T value);

    //! Set or get an option in string representation
    std::string& operator[](const std::string& name);


    //! Delete an option
    void delete_option(const std::string& name);

    //! Clear all options
    void clear(void);

    //! operator to add options
    ModelOptions& operator+=(const ModelOptions& rhs);

    //! operator to add options
    ModelOptions& operator+=(const std::map<const std::string,
        std::string>& rhs);

    //! Print all options for debugging
    void print_all(void) const;


  private:

    //! typedef for convenience
    typedef std::map<const std::string, std::string> OptionsMap;
    
    //! The map holding all options
    OptionsMap _options;

};




//
// inline methods
//


inline
void
ModelOptions::delete_option(const std::string& name)
{
  _options.erase(name);
}



inline
std::string&
ModelOptions::operator[](const std::string& name)
{
  return _options[name];
}



inline
bool
ModelOptions::find_option(const std::string& name) const
{
  bool res = true;

  OptionsMap::const_iterator it(_options.find(name));
  
  if (it == _options.end())
    res = false;

  return res;
}



inline
void
ModelOptions::clear(void)
{
  _options.clear();
}



template <typename T>
inline
T
ModelOptions::get_option(const std::string& name, T default_value) const
{
  OptionsMap::const_iterator it = _options.find(name);

  if (it != _options.end())
    return Utils::convert<T>(it->second);
  else
    return default_value;
}


inline
bool
ModelOptions::is_empty(void)
{
  return _options.empty();
}


#endif // _MODELOPTIONS_H_
