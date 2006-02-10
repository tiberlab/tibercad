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

    //! Get the identifier of the property object
    /*!
     * The identifier is a unique string which identifies each type
     * of physical properties.
     */
    virtual const std::string get_id(void) const = 0;

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


  protected:
    
    //! The empty constructor
    /*!
     * \c PhysicalProperties should not be instantiated directly
     */
    PhysicalProperties(void);

    //! The empty destructor
    /*!
     * should be implemented in the derived classes if needed
     */
    virtual ~PhysicalProperties(void) {};

  private:

    //! The material this properties belong to
    const Material* _material;
};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
PhysicalProperties::PhysicalProperties(void)
  : _material(0)
{
}

inline void
PhysicalProperties::set_material(const Material* material)
{
  _material = material;
}


#endif // _PHYSICALPROPERTIES_H_
