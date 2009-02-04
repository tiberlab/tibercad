// $Id$

#ifndef _DATABASE_H_
#define _DATABASE_H_

#include "getpot.h"
#include <string>

class Database
{

  public:

    //! Constructor
    Database(void);

    //! Set the search path for the material files
    void set_search_path(const std::string& path);

    //! Get the search path for material files
    const std::string& get_search_path(void) const;

    //! Set the system wide default search path
    static void set_default_search_path(const std::string& path);

    //! Set the material
    /*!
     * \param material the name of the material
     * \param datafile the file containing the material parameters
     *
     * If \c datafile is empty, the material file is looked for in the
     * search path.
     */
    void set_material(const std::string& material,
        const std::string& datafile = "");

    //! Get the material
    const std::string& get_material(void) const;


    //! \deprecated
    const std::string get_data_file(void) const
      { return _datafile; };

    //! Set the name of the model
    void set_section(const std::string& section);

    //! Get the name of the model
    const std::string& get_section(void) const;

    //! Returns \c true if the material is an alloy
    bool is_alloy(const std::string& material) const;

    //! Returns the alloy components
    void get_alloy_components(const std::string& alloy,
        std::string& comp_A, std::string& comp_B);

    //! Returns the alloy components
    void get_alloy_components(std::string& comp_A,
        std::string& comp_B) const;


    //! Checks if a given variable is present in the databas
    bool has_variable(const std::string& variable) const;


    //! Get data of some type
    template <typename T>
    T get(const std::string& variable, T default_value,
        bool required = false) const;


    //! Get string data
    std::string get(const std::string& variable,
        const std::string& default_value, bool required = false) const;

    //! Get string data
    std::string get(const std::string& variable,
        const char* default_value, bool required = false) const;


    //! Get data array/vector
    /*!
     * If \c variable is not present in the database, the input vector
     * will be left unchanged. If a non-empty vector is provided, it will
     * throw an exception if it cannot find the same number of elements in the
     * database.
     */
    template <typename T>
    void get(const std::string& variable, std::vector<T>& data,
        bool required = false) const;

    //! Get data matrix
    /*!
     * If \c variable is not present in the database, the input matrix
     * will be left unchanged. If a non-empty matrix is provided, it will
     * throw an exception if the database does not provide the same data
     * structure.
     */
    template <typename T>
    void get(const std::string& variable,
        std::vector<std::vector<T> >& data, bool required = false) const;



  private:

    //! The search path for the material files
    std::string _path;

    //! The default search path
    static std::string _default_path;

    //! The name of the material currently processed
    std::string _material;

    //! The name of the model currently processed
    std::string _section;

    //! The currently used data file
    std::string _datafile;

    //! The file parser
    GetPot _file;

    //! Check the data file
    bool check_data_file(const std::string& name) const;

    //! Get the name of the currently used data file
    const std::string get_data_file(const std::string& material) const;

    //! Check for a variable and throw exception if it is not found
    void require_variable(const std::string& variable) const;

};


//
// inline methods
//

inline
Database::Database(void)
  : _path(""),
    _material(""),
    _section(""),
    _datafile("")
{
}

inline
const std::string&
Database::get_material(void) const
{
  return _material;
}


inline
const std::string&
Database::get_section(void) const
{
  return _section;
}


inline
bool
Database::has_variable(const std::string& variable) const
{
  return _file.have_variable(variable.c_str());
}



inline
std::string
Database::get(const std::string& variable,
    const std::string& default_value, bool required) const
{
  if (required) require_variable(variable);
  return _file(variable.c_str(), default_value);
}


inline
std::string
Database::get(const std::string& variable,
    const char* default_value, bool required) const
{
  if (required) require_variable(variable);
  return _file(variable.c_str(), std::string(default_value));
}



template <typename T>
inline
T
Database::get(const std::string& variable, T default_value,
    bool required) const
{
  if (required) require_variable(variable);
  return _file(variable.c_str(), default_value);
}



#endif // _DATABASE_H_
