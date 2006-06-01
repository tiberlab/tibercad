#ifndef _BOUNDARYPROPERTIES_H_
#define _BOUNDARYPROPERTIES_H_

#include <string>

#include "types.h"


//! Boundary properties for a certain type of simulation.
/*!
 *  This is the base class for boundary properties. Every solver
 *  module should implement its own class derived from this one to hold
 *  the properties needed for calculations.
 */
class BoundaryProperties
{

  public:

    //! Get the name of the boundary
    /*!
     * The name is a string which identifies a boundary.
     */
    const std::string get_name(void);

    //! Get ID of the boundary
    /*!
     * The ID of a boundary should be a unique number
     */
    ID get_id(void) const;

    //! Set the name of the boundary
    /*!
     * Set the logical name of the boundary
     */
    void set_name(const std::string& name);


  protected:
    
    //! The empty constructor
    /*!
     * \c BoundaryType should not be instantiated directly
     *
     * \param id the unique ID for the boundary
     * \param name the name of the boundary
     */
    BoundaryProperties(ID id, const std::string name = "");

    //! The empty destructor
    /*!
     * should be implemented in the derived classes if needed
     */
    virtual ~BoundaryProperties(void) {};

  private:

    //! The identifier of a boundary object
    /*!
     * The identifier is assumed to be uniqe for every instance of
     * \c Boundary
     */
    Id _id;

    //! The name of a boundary object
    /*!
     * This is a logical name for the boundary as e.g. "gate"
     */
    const std::string _name;

    //! The list of element sides which lie on this boundary
    /*!
     * The list contains all element sides of the level 0 mesh
     * which are lying on this boundary
     */
    std::set<ElemSide> _boundary_elements;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
BoundaryProperties::Boundary(ID id, const std::string name = "")
  : _id(id), _name(name)
{
}


inline
BoundaryProperties::set_name(const std::string name)
{
  _name = name;
}


#endif // _BOUNDARYPROPERTIES_H_
