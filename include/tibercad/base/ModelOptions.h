/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file ModelOptions.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_MODELOPTIONS_H
#define TC_MODELOPTIONS_H

#include "tibercad/utils/Utils.h"
#include "tibercad/base/tiber_dll.h"

#include <map>
#include <set>
#include <vector>
#include <string>

namespace libMesh
{
  class Point;
}


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
    ModelOptions(void);


    //! The copy constructor
    ModelOptions(const ModelOptions& other);


    //! Constructor which takes a map as argument
    ModelOptions(std::map<const std::string, std::string> options);


    //! The destructor
    ~ModelOptions(void) {};


    //! Check if it is empty
    bool is_empty(void) const;


    //! Set the input file block keyword
    void set_key(const std::string& key);


    //! Get the input file block keyword
    const std::string& get_key(void) const;


    //! Set the name
    /*!
     * Sets also the option \c name to the given value
     */
    void set_name(const std::string& modifier);


    //! Get the name
    /*!
     * The name is set from the second keyword in the input file
     */
    const std::string& get_name(void) const;

    //! Get the index of the model
    int get_index(void) const;


    //! Get an option as string as given in input file
    std::string get_raw_option_string(const std::string& name) const;


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


    //! Get an option which is a set of values (of the same type)
    /*!
     * \param name the name of the option
     * \param data the set, where the values will be stored. \c data can
     * contain default values.
     */
    template <typename T>
    void get_option(const std::string& name, std::set<T>& data) const;


    //! Get an option which is a three component vector
    /*!
     * \param name the name of the option
     * \param vec the vector, where the values will be stored.
     */
    void get_option(const std::string& name, libMesh::RealVectorValue& vec) const;


    //! Get an option which is a point
    /*!
     * \param name the name of the option
     * \param point the point, where the values will be stored.
     */
    void get_option(const std::string& name, libMesh::Point& vec) const;

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
     * \param value the value to set/home/maufder/Work/TiberCAD/TransferToGit/include/base/common/ModelOptions.h
     */
    template <typename T>
    void set_option(const std::string& name, const T& value);


    //! Set an option providing a const char*
    void set_option(const std::string& name, const char* value);


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


    //! Assignment operator
    ModelOptions& operator=(const ModelOptions& rhs);


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
     * \note The model name doesn't have to be unique.
     */
    void add_submodel(const std::string& name, const ModelOptions& options);


    //! Delete a certain submodel
    /*!
     * \param it the iterator to the submodel options to be deleted
     */
    void delete_submodel(submodel_iterator it);


    //! Delete all submodels of a certain kind
    /*!
     * \param name the name of the model to be deleted
     */
    void delete_submodels(const std::string& name);


    //! Delete all submodels
    void delete_all_submodels(void);


    //! Check if there is a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return \c true if found, \c false otherwise
     */
    bool has_submodel(const std::string& name) const;


    //! Get the iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the const iterator for the first appearance of the model
     */
    submodel_iterator submodels_begin(const std::string& name);


    //! Get the past-the-end iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the past-the-end iterator for the model
     */
    submodel_iterator submodels_end(const std::string& name);


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


    //! Get the const iterator for the first submodel
    const_submodel_iterator submodels_begin(void) const;


    //! Get the past-the-end iterator for the submodels
    const_submodel_iterator submodels_end(void) const;


    //! Check if there are unused options
    /*!
     * mode = 0, do nothing
     * mode = 1, warn
     * mode = 2, throw exception
     *
     * return \c true if unused options are found
     */
    bool check_unused(int mode = 1) const;


  private:

    //! The block keyword
    /*!
     * This is the first keyword in the input file
     */
    std::string _key;

    //! A name
    /*!
     * This is set from the second keyword in the input file
     */
    std::string _name;

    //! A counter used to count blocks of the same type
    /*!
     * This can be used to identify submodels based on the order they are defined
     */
    int _index;

    //! The map holding all options
    OptionsMap _options;

    //! A map containing submodels
    SubmodelMap _submodels;

    //! Get the iterator for an option
    OptionsMap::const_iterator _find(const std::string& name) const;


    //mutable std::list<OptionsMap::const_iterator> _used;
    mutable std::set<std::string> _used;

};




//
// inline methods
//

inline
ModelOptions::ModelOptions(void) :
  _key(""),
  _name(""),
  _index(-1)
{

}





inline
const std::string&
ModelOptions::get_name(void) const
{
  return _name;
}


inline
void
ModelOptions::set_key(const std::string& key)
{
  _key = key;
}


inline
const std::string&
ModelOptions::get_key(void) const
{
  return _key;
}


inline
int
ModelOptions::get_index(void) const
{
  return(_index);
}


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

  if (_find(name) == _options.end())
    res = false;

  return res;
}


inline
void
ModelOptions::set_option(const std::string& name, const char* value)
{
  _options[name] = value;
}



template <typename T>
inline
T
ModelOptions::get_option(const std::string& name, T default_value) const
{
  OptionsMap::const_iterator it(_find(name));

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
ModelOptions::submodel_iterator
ModelOptions::submodels_begin(const std::string& name)
{
  return _submodels.lower_bound(name);
}



inline
ModelOptions::submodel_iterator
ModelOptions::submodels_end(const std::string& name)
{
  return _submodels.upper_bound(name);
}






inline
bool
ModelOptions::has_submodel(const std::string& name) const
{
  bool ans = true;
  if (_submodels.find(name) == _submodels.end())
    ans = false;

  return ans;
}


/*
inline
ModelOptions::const_option_iterator
ModelOptions::options_begin(void) const
{
  return _options.begin();
}


inline
ModelOptions::const_option_iterator
ModelOptions::options_end(void) const
{
  return _options.end();
}
*/

#endif // TC_MODELOPTIONS_H
