// $Id$

#ifndef _DATABASE_H_
#define _DATABASE_H_

#include <string>
#include <vector>

class GetPot;

class Database
{

  public:

    //! A type to define different alloy mixing methods
    enum AlloyMixing
    {
      NONE = 0, //!< do not do any mixing at all
      VCA       //!< do virtual crystal approximation
    };

    //! Default constructor
    Database(void);

    //! Constructor taking arguments
    Database(const std::string& material,
        const std::string& datafile = "");

    //! Copy constructor
    Database(const Database& other);

    //! Destructor
    ~Database(void);

    //! Assignement operator
    Database& operator=(const Database& rhs);


    //! Set the search path for the material files
    static void set_search_path(const std::string& path);

    //! Get the search path for material files
    static const std::string& get_search_path(void);

    //! Set the system wide default search path
    static void set_default_search_path(const std::string& path);

    //! Get the system wide default search path
    static std::string& get_default_search_path(void);

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


    //! Set the data file
    void set_data_file(const std::string& file);


    //! Get the data file name
    const std::string get_data_file(void) const;



    //! Call this after you have finished reading
    void close(void);


    //! Set the alloy mixing type
    void set_alloy_mixing(AlloyMixing type);


    //! Set the alloy mixing type
    AlloyMixing get_alloy_mixing(void) const;


    //! Set the name of the model
    /*!
     * Call this method always before reading data.
     */
    void set_section(const std::string& section);

    //! Get the name of the model
    const std::string& get_section(void) const;

    //! Returns \c true if the material is an alloy
    bool is_alloy(void) const;


    //! Returns the alloy components
    void get_alloy_components(std::string& comp_A,
        std::string& comp_B) const;


    //! Returns the alloy components
    void get_alloy_components(std::vector<std::string>& comp) const;


    //! Get the number of components
    size_t get_number_of_components(void) const;


    //! Get the database of the i-th component
    /*!
     * \param i the component index, beginning from 1
     */
    const Database& get_component_database(size_t i) const;


    //! For an alloy, set the composition (fraction for each component)
    void set_alloy_composition(std::vector<double>& fractions);


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
    static std::string _path;

    //! The default search path
    static std::string _default_path;

    //! The name of the material
    std::string _material;

    //! The name of the model currently processed
    mutable std::string _section;

    //! The currently used data file
    std::string _datafile;

    //! The file parser
    /*!
     * The actual file pointer does not make part of the
     * internal state of the database, and we want to assign it
     * on the fly when accessing the database.
     */
    mutable GetPot* _file;


    //! \c true if this is the database of an alloy
    bool _is_alloy;


    //! The type of alloy mixing
    AlloyMixing _mixing_type;


    //! The databases of the constituents for an alloy
    std::vector<Database> _comp_db;


    //! The component fractions
    std::vector<double> _comp_fractions;


    //! Check the data file
    bool check_data_file(const std::string& name) const;

    //! Get the automatic data file name
    const std::string get_data_file(const std::string& material) const;

    //! Open the database
    /*!
     * This method is called by all other methods that
     * access the database
     */
    void open(void) const;

    //! Does the real opening of the databas
    void do_open(void) const;

    //! Check for a variable and throw exception if it is not found
    void require_variable(const std::string& variable) const;

};


//
// inline methods
//

inline
Database::Database(void)
  : _material(""),
    _section(""),
    _datafile(""),
    _file(NULL),
    _is_alloy(false),
    _mixing_type(VCA)
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
std::string&
Database::get_default_search_path(void)
{
  return _default_path;
}


inline
void
Database::set_data_file(const std::string& file)
{
  _datafile = file;
}

inline
const
std::string
Database::get_data_file(void) const
{
  return _datafile;
}


inline
void
Database::open(void) const
{
  if (_file == NULL)
    do_open();
}


inline
void
Database::set_alloy_mixing(AlloyMixing mixing_type)
{
  _mixing_type = mixing_type;
}


inline
Database::AlloyMixing
Database::get_alloy_mixing(void) const
{
  return _mixing_type;
}


inline
bool
Database::is_alloy(void) const
{
  return _is_alloy;
}


inline
size_t
Database::get_number_of_components(void) const
{
  return _comp_db.size();
}


inline
const Database&
Database::get_component_database(size_t i) const
{
  return _comp_db.at(i);
}




#endif // _DATABASE_H_
