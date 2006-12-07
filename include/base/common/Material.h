// $Id$

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "PhysicalModel.h"

#include "point.h"

#include <string>
#include <map>
#include <vector>

// forward declarations
class Database;
class RotatedCrystal;

//! Contains all needed data for a material
/*!
 * For any material/structure combination that is used for a simulation
 * a \c Material object is built which contains all needed physical
 * model
 */
class Material
{

  public:

    //! Destructor
    /*!
     * Deletes all \c PhysicalModel objects
     */
    virtual ~Material(void);

    //! Set the database to be used
    /*!
     * The database could in principle be overriden from the options
     */
    static void set_database(Database& database);

    //! Create a material with name \c name
    static Material* create(const std::string& name);

    //! Create a material with name \c name and options
    static Material* create(const std::string& name,
        const ModelOptions& options);

    //! Initialize the material
    /*!
     * Read all needed material data from the database
     * It calls the \c read_database() of each \c PhysicalModel
     * object
     */
    void init(void);

    //! Add new physical model
    /*!
     * Add a new \c PhysicalModelInterface object for this material
     *
     * \param model the model to add
     * \param simulator_id the id of the simulator this model is used for
     */
    void add_model(PhysicalModel* model, ID simulator_id);

    //! Get physical model for simulator with ID id
    /*!
     *
     * It will return the \c NULL pointer if the requested model
     * is not in the list.
     * \c id is a unique ID assigned to each simulator.
     *
     * \param id the identifier for simulator
     * \return a pointer to the model object
     */
    PhysicalModel* get_model(ID id) const;

    //! Get the material name
    const std::string& get_name(void) const;

    //! Get the crystal structure
    const std::string& get_structure(void) const;

    //! Get a reference to the RotatedCrystal
    const RotatedCrystal& get_rotated_crystal(void) const;

    //! Get a reference to the database
    const Database& get_database(void) const;

  protected:

    //! a typedef for convenience
    typedef std::map<ID, PhysicalModel*> ModelMap;

    
    //! Construct a material with a given structure
    /*!
     * At construction one has to specify the material name
     *
     * \param name the name of the Material
     */
    Material(const std::string& name);

    //! The real init function
    /*!
     * This one gets called from init(const Database& db)
     */
    virtual void do_init(void);

    //! Set the model options
    void set_options(const ModelOptions& options);

    //! Set the structure
    void set_structure(const std::string& structure);

    //! Get the options
    ModelOptions& get_options(void);
    
    //! Get a writable reference to the database
    Database& get_database(void);

    //! Get an iterator to the first model
    ModelMap::iterator models_begin(void);

    //! Get an iterator to the last model
    ModelMap::iterator models_end(void);


  private:

    //! The material name
    /*!
     * The name of the material as GaAs, Si etc.
     */
    const std::string _name;

    //! The crystal structure
    /*!
     * The crystal structure as wz, zb etc
     */
    std::string _structure;

    //! The RotatedCrystal object which
    RotatedCrystal* _rotated_crystal;

    //! The map containing all \c PhysicalModelInterface objects
    /*!
     * This map containes the physical model of any simulation type
     * requested.
     */
    ModelMap _models;

    //! The default database to be used
    static Database* _database;

    //! Options for this material
    ModelOptions _options;


};


//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
Material::Material(const std::string& name)
: _name(name), _structure("zb")
{
}


inline
PhysicalModel*
Material::get_model(ID id) const
{
  const ModelMap::const_iterator end(_models.end());
  ModelMap::const_iterator it(_models.find(id));
  if (it != end)
    return it->second;
  else
    return NULL;
}


inline
const std::string&
Material::get_name(void) const
{
  return  _name;
}


inline
const std::string&
Material::get_structure(void) const
{
  return  _structure;
}


inline
const RotatedCrystal&
Material::get_rotated_crystal(void) const
{
  return *_rotated_crystal;
}


inline
void
Material::set_structure(const std::string& structure)
{
  _structure = structure;
}


inline
const Database&
Material::get_database(void) const
{
  return *_database;
}


inline
Database&
Material::get_database(void)
{
  return *_database;
}


inline
void
Material::init(void)
{
  assert(_database != NULL);
  
  do_init();
}


inline
void
Material::set_database(Database& database)
{
  _database = &database;
}

inline
void
Material::set_options(const ModelOptions& options)
{
  _options = options;
}

inline
ModelOptions&
Material::get_options(void)
{
  return _options;
}


inline
Material::ModelMap::iterator
Material::models_begin(void)
{
  return _models.begin();
}


inline
Material::ModelMap::iterator
Material::models_end(void)
{
  return _models.end();
}


#endif // _MATERIAL_H_
