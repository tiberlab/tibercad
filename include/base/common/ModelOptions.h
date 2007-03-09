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

    //! typedef for convenience
    typedef std::map<const std::string, std::string> OptionsMap;
    
    //! typedef for the map of submodels
    typedef std::multimap<const std::string, ModelOptions> SubmodelMap;

    //! An iterator to iterate over the submodels
    typedef SubmodelMap::iterator submodel_iterator;

    /*! \copydoc submodel_iterator */
    typedef SubmodelMap::const_iterator const_submodel_iterator;


    //! The default constructor
    ModelOptions(void) {};

    
    //! Constructor which takes a map as argument
    ModelOptions(std::map<const std::string, std::string> options);

    
    //! The destructor
    ~ModelOptions(void) {};

    
    //! Check if it is empty
    bool is_empty(void) const;

    
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
     * \param name the name of the option
     * \param value the value to set
     */
    template <typename T>
    void set_option(const std::string& name, const T& value);


    //! Set an vector option
    /*!
     * \param name the name of the option
     * \param value the vector with the values to set
     */
    template <typename T>
    void set_option(const std::string& name, const std::vector<T>& value);
    
    
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


    //! Add a submodel
    /*!
     * \param name the name of the model to be added
     * \param options the ModelOptions
     *
     * \note { The model name doesn't have to be unique. }
     */
    void add_submodel(const std::string& name, const ModelOptions& options);


    //! Check if there is a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return \c true if found, \c false otherwise
     */
    bool has_submodel(const std::string& name);


    //! Get the const iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the const iterator for the first appearance of the model
     */
    const_submodel_iterator submodels_begin(const std::string& name) const;


    //! Get the past-the-end iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the past-the-end iterator for the model
     */
    const_submodel_iterator submodels_end(const std::string& name) const;


    //! Get the iterator for the first submodel
    submodel_iterator submodels_begin(void);

    
    //! Get the past-the-end iterator for the submodels
    submodel_iterator submodels_end(void);

    
    //! Get the iterator for the first submodel
    const_submodel_iterator submodels_begin(void) const;

    
    //! Get the past-the-end iterator for the submodels
    const_submodel_iterator submodels_end(void) const;




  private:
    
    //! The map holding all options
    OptionsMap _options;

    //! A map containing submodels
    SubmodelMap _submodels;

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
ModelOptions::is_empty(void) const
{
  return _options.empty();
}



inline
ModelOptions::const_submodel_iterator
ModelOptions::submodels_begin(void) const
{
  return _submodels.begin();
}



inline
ModelOptions::const_submodel_iterator
ModelOptions::submodels_end(void) const
{
  return _submodels.end();
}




inline
ModelOptions::submodel_iterator
ModelOptions::submodels_begin(void)
{
  return _submodels.begin();
}




inline
ModelOptions::submodel_iterator
ModelOptions::submodels_end(void)
{
  return _submodels.end();
}



inline
ModelOptions::const_submodel_iterator
ModelOptions::submodels_begin(const std::string& name) const
{
  return _submodels.lower_bound(name);
}



inline
ModelOptions::const_submodel_iterator
ModelOptions::submodels_end(const std::string& name) const
{
  return _submodels.upper_bound(name);
}




inline
void
ModelOptions::add_submodel(const std::string& name,
    const ModelOptions& options)
{
  _submodels.insert(SubmodelMap::value_type(name, options));
}



inline
bool
ModelOptions::has_submodel(const std::string& name)
{
  bool ans = true;
  if (_submodels.find(name) == _submodels.end())
    ans = false;

  return ans;
}



    
#endif // _MODELOPTIONS_H_
