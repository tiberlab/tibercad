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

    //! Get the name of the currently used data file
    const std::string get_data_file(const std::string& material) const;


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


    //! Get double data
    double get(const std::string& variable, double default_value);

    //! Get integer data
    int get(const std::string& variable, int default_value);

    //! Get bool data
    bool get(const std::string& variable, bool default_value);

    //! Get string data
    std::string get(const std::string& variable,
        const std::string& default_value);

    //! Get string data
    std::string get(const std::string& variable,
        const char* default_value);



  private:

    //! The search path for the material files
    std::string _path;

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
double
Database::get(const std::string& variable, double default_value)
{
  return _file(variable.c_str(), default_value);
}


inline
int
Database::get(const std::string& variable, int default_value)
{
  return _file(variable.c_str(), default_value);
}


inline
bool
Database::get(const std::string& variable, bool default_value)
{
  return _file(variable.c_str(), default_value);
}


inline
std::string
Database::get(const std::string& variable, const std::string& default_value) 
{
  return _file(variable.c_str(), default_value);
}


inline
std::string
Database::get(const std::string& variable, const char* default_value) 
{
  return _file(variable.c_str(), std::string(default_value));
}


#endif // _DATABASE_H_
