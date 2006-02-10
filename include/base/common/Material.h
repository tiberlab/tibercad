// $Id$

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <string>
#include <map>

#include "PhysicalProperties.h"

// forward declarations
class Dummy;

//! Contains all needed data for a material
/*!
 * For any material/structure combination that is used for a simulation
 * a \c Material object is built which contains all needed physical
 * properties
 */
class Material
{

  public:

    //! Construct a material with a given structure
    /*!
     * At construction one has to specify at least the material name and
     * optionally also the structure. The default for the structure is 
     * zincblende.
     *
     * \param name the name of the Material
     * \param structure the crystal structure
     */
    Material(const std::string& name, const std::string& structure = "zb");

    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    ~Material(void);

    //! Initialize the material
    /*!
     * Read all needed material data from the database
     * It calls the \c read_database() of each \c PhysicalProperties
     * object
     *
     * \param database the database to read from
     */
    void init(const Dummy& database);

    //! Add new physical properties
    /*!
     * Add a new \c PhysicalProperties object for this material
     *
     * \param properties the properties to add
     */
    void add_properties(PhysicalProperties* properties);

    //! Get physical properties of a given type
    /*!
     * Get the physical properties of type \c id.
     * It will return the \c NULL pointer if the requested properties
     * are not in the list.
     *
     * \param id the identifier for the set of properties
     * \return a const pointer to the property object
     */
    const PhysicalProperties* get_properties(const std::string& id) const;


  private:

    //! a typedef for convenience
    typedef std::map<const std::string, PhysicalProperties*> PropertyMap;

    //! the material name
    /*!
     * The name of the material as GaAs, Si etc.
     */
    const std::string _name;

    //! the crystal structure
    /*!
     * The crystal structure as wz, zb etc
     */
    const std::string _structure;

    //! the map containing all \c PhysicalProperties objects
    /*!
     * This map containes the physical properties of any simulation type
     * requested.
     */
    PropertyMap _properties;

};


//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
Material::Material(const std::string& name, const std::string& structure)
  : _name(name), _structure(structure)
{
}

inline
Material::~Material(void)
{
}

inline const PhysicalProperties*
Material::get_properties(const std::string& id) const
{
  const PropertyMap::const_iterator end = _properties.end();
  PropertyMap::const_iterator it = _properties.find(id);
  if (it != end)
    return it->second;
  else
    return NULL;
}

inline void
Material::add_properties(PhysicalProperties* properties)
{
  const std::string& id = properties->get_id();
  _properties[id] = properties;
}

inline void
Material::init(const Dummy& database)
{
  PropertyMap::iterator it = _properties.begin();
  const PropertyMap::const_iterator end = _properties.end();

  for ( ; it != end; ++it)
  {
    (it->second)->set_material(this);
    (it->second)->read_database(database);
  }
}



#endif // _MATERIAL_H_
