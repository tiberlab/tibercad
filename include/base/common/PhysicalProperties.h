#ifndef _PHYSICALPROPERTIES_H_
#define _PHYSICALPROPERTIES_H_

#include "TypeDefs.h"

#include <string>
#include <map>
#include <iostream>

// forward declarations
class Dummy;
class Material;



//! Some physical properties of a material.
/*!
 *  This is the base class for physical material properties. Every solver
 *  module should implement its own class derived from this one to hold
 *  the properties needed for calculations.
 */
class PhysicalProperties
{

 public:


  //! The empty destructor
  /*!
   * should be implemented in the derived classes if needed
   */
  virtual ~PhysicalProperties(void) {};

  //! Get the name of the property object
  /*!
   * The name is a unique string which identifies each type
   * of physical properties.
   */
  const std::string& get_name(void) const;


  //! Get the identifier of the property object
  /*!
   * The identifier is a unique numerical ID which identifies each type
   * of physical properties.
   */
  ID get_id(void) const;

  //! Get the identifier of the property object \c name
  /*!
   * The identifier is a unique numerical ID which identifies each type
   * of physical properties.
   *
   * \param name the name of the property object
   */
  static ID get_id(const std::string& name);



  //! Set a pointer to the \c Material object this properties belong to
  /*!
   * \param material a constant pointer to a material
   */
  void set_material(const Material* material);

  //! Read the properties from the database
  /*!
   * Reads all needed physical properties from the database.
   * \param db a const reference to the database
   */
  virtual void read_database(const Dummy& db) = 0;  //  for  compatibility with  svn !!
  //  virtual void read_database(DataBaseCall& db) = 0;



  //! Read the bowing parameters  from the database
  /*!
   * Read the bowing parameters for an alloy material from the database.
   * \param db a const reference to the database
   */
  virtual void read_database_bowing_parameters(const Dummy& db ) {};

   
  //! Constructs  \c PhysicalProperties for  an  alloy.
  /*!
   * Constructs  \c PhysicalProperties for  an  alloy, basing on
   * properties of  component materials
   * and on alloy composition \param molar_fraction.
   */ 
    
  virtual void set_properties_alloy(const  PhysicalProperties* prop_comp1, 
				    const  PhysicalProperties* prop_comp2, 
				    double molar_fraction) {};





  //! Returns pointer to the  material these properties belong  to.
  /*!
   * 
   */
  const Material*  get_material(void) const;

 protected:


  //! The empty constructor
  /*!
   * \c PhysicalProperties should not be instantiated directly
   *
   * \param id a unique string which identifies the type of properties
   */
  PhysicalProperties(const std::string& name);





 private:

  //! For convenience
  typedef std::map<const std::string, ID> PropertyMap;

  //! The name of this properties object
  const std::string _name;

  //! The id of this properties object
  ID _id;

  //! The material this properties belong to
  const Material* _material;

  //! A map which assigns numerical IDs to a property name
  static PropertyMap _property_map;

  //! The ID to be given to the next property type
  static ID _max_id;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
PhysicalProperties::PhysicalProperties(const std::string& name)
  : _name(name), _id(0), _material(0)
{
  if (_property_map.find(name) == _property_map.end())
  {
    _id = _max_id;
    _property_map[name] = _id;
    ++_max_id;
  }
  //std::cerr << name << " " << get_id(name) << "\n";
}

inline
void
PhysicalProperties::set_material(const Material* material)
{
  _material = material;
}

inline
const std::string&
PhysicalProperties::get_name(void) const
{
  return _name;
}

inline
ID
PhysicalProperties::get_id(const std::string& name)
{
  ID id = 0;
  std::map<const std::string, ID>::iterator it(_property_map.find(name));
  if (it != _property_map.end())
    id = it->second;

  return id;
}

inline
ID
PhysicalProperties::get_id(void) const
{
  return _id;
}

inline
const Material*
PhysicalProperties::get_material(void) const
{
  return  _material;
}





#endif // _PHYSICALPROPERTIES_H_


