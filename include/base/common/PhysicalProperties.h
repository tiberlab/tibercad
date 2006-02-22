#ifndef _PHYSICALPROPERTIES_H_
#define _PHYSICALPROPERTIES_H_

#include <string>

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
    
    //! Get the identifier of the property object
    /*!
     * The identifier is a unique string which identifies each type
     * of physical properties.
     */
    const std::string& get_id(void) const;


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
  virtual void read_database(const Dummy& db) = 0;
    
    
  virtual void read_database_bowing_parameters(const Dummy& db ) {};
    
    
  virtual void set_properties_alloy(const  PhysicalProperties* prop_comp1, 
				    const  PhysicalProperties* prop_comp2, 
				    double molar_fraction) {};


 protected:
    

    //! The empty constructor
    /*!
     * \c PhysicalProperties should not be instantiated directly
     */
    PhysicalProperties(const std::string& id);





 private:


    //! The id of this properties
    const std::string _id;

    //! The material this properties belong to
    const Material* _material;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
PhysicalProperties::PhysicalProperties(const std::string& id)
  : _id(id), _material(0)
{
}

inline
void
PhysicalProperties::set_material(const Material* material)
{
  _material = material;
}

inline
const std::string&
PhysicalProperties::get_id(void) const
{
  return _id;
}


#endif // _PHYSICALPROPERTIES_H_
