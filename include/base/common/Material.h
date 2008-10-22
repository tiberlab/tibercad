// $Id$

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "PhysicalModel.h"

// LibMesh includes
#include "point.h"

// C++ includes
#include <string>
#include <map>
#include <set>
#include <vector>

// forward declarations
class Dopant;
class Database;
class RotatedCrystal;


//! Contains all needed data for a material
/*!
 * For any material/structure combination that is used for a simulation
 * a \c Material object is built which contains all needed physical
 * models. Every simulation will have exactly one model in the models list
 * which it can use for its calculations. Additionally, a Material object
 * contains also the list of donors and acceptors (see Dopant) and a
 * RotatedCrystal object.
 */
class Material
{

  public:

    //! An iterator to iterate over all dopants
    typedef std::set<Dopant*>::iterator dopant_iterator;

    //! A const iterator to iterate over all dopants
    //typedef std::set<Dopant*>::const_iterator const_dopant_iterator;


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


    //! Tells if this material is an alloy
    bool is_alloy(void) const;

    
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

    
    //! Add a dopant
    void add_dopant(Dopant* dopant);

    
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
    
    
    //! Get a writable reference to the database
    Database& get_database(void);

    
    //! Get the options
    ModelOptions& get_options(void);

    
    //! Get the options
    const ModelOptions& get_options(void) const;

    
    //! Get the total n-doping
    double get_total_donor_density(void) const;

    
    //! Get the total p-doping
    double get_total_acceptor_density(void) const;


    //! Get the total doping density
    /*!
     * The return value is \f$N_d + N_a\f$
     */
    double get_total_doping_density(void) const;


    //! Get the total net doping density
    /*!
     * The return value is \f$N_d - N_a\f$
     */
    double get_net_doping_density(void) const;


    //! Get the first iterator for the donors
    dopant_iterator donors_begin(void) const;

    //! Get the past-the-end iterator for the donors
    dopant_iterator donors_end(void) const;


    //! Get the first iterator for the acceptors
    dopant_iterator acceptors_begin(void) const;

    //! Get the past-the-end iterator for the acceptors
    dopant_iterator acceptors_end(void) const;

    
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

    
    //! Get a writable reference to the RotatedCrystal
    RotatedCrystal& get_crystal(void);

    
    //! Get an iterator to the first model
    ModelMap::iterator models_begin(void);

    
    //! Get an iterator to the last model
    ModelMap::iterator models_end(void);
    

    //! True if this is an alloy
    bool _is_alloy;

    
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

    
    //! The RotatedCrystal object
    RotatedCrystal* _rotated_crystal;


    //! The map containing all \c PhysicalModelInterface objects
    /*!
     * This map containes the physical model of any simulation type
     * requested.
     */
    ModelMap _models;


    //! The list of donors
    std::set<Dopant*> _donors;

    //! The list acceptors
    std::set<Dopant*> _acceptors;

    
    //! The default database to be used
    static Database* _database;

    
    //! Options for this material
    ModelOptions _options;


    //! Clear all doping
    void clear_doping(void);


    //! A flag to tell if the material is already initialized
    bool _is_initialized;


};


//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
Material::Material(const std::string& name)
  : _is_alloy(false),
    _name(name),
    _structure("zb"), 
    _rotated_crystal(NULL),
    _is_initialized(false)
{
}

inline
bool
Material::is_alloy(void) const
{
  return _is_alloy;
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
RotatedCrystal&
Material::get_crystal(void)
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
  
  if (!_is_initialized)
  {
    do_init();
    _is_initialized = true;
  }
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
const ModelOptions&
Material::get_options(void) const
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


inline
double
Material::get_total_doping_density(void) const
{
  return (get_total_donor_density() + get_total_acceptor_density());
}


inline
double
Material::get_net_doping_density(void) const
{
  return (get_total_donor_density() - get_total_acceptor_density());
}


inline
Material::dopant_iterator
Material::donors_begin(void) const
{
  return _donors.begin();
}


inline
Material::dopant_iterator
Material::donors_end(void) const
{
  return _donors.end();
}


inline
Material::dopant_iterator
Material::acceptors_begin(void) const
{
  return _acceptors.begin();
}


inline
Material::dopant_iterator
Material::acceptors_end(void) const
{
  return _acceptors.end();
}





#endif // _MATERIAL_H_
