// $Id$


#ifndef _SOLUTIONDESCRIPTOR_H_
#define _SOLUTIONDESCRIPTOR_H_


//! A structure describing the properties of a solution
struct SolutionDescriptor
{

  public:
    enum Type
    {
      REAL,     //!< a real value
      COMPLEX,  //!< a complex value, ordered as (real,imag)
      VECTOR,   //!< a real vector with three components (x, y, z)
      TENSOR    //!< a real symmetric tensor of second rank ( TODO )
    };


    enum Location
    {
        NODAL,
        CELL,
        ATOM,
        GLOBAL
    };

    //! Default constructor
    /*!
     * Sets the ID to \c INVALID_ID
     */
    SolutionDescriptor(void) : _id(INVALID_ID) {}

    SolutionDescriptor(const std::string& name, ID id,
        Type type, Location location, const std::string& units = "");


    //! Get the ID
    ID id(void) const { return _id; }

    //! Get the name
    const std::string& name(void) const { return _name; }

    //! Get the type
    Type type(void) const { return _type; }

    //! Get the location
    Location location(void) const { return _location; }

    //! Get the units
    const std::string& units(void) const { return _units; }


  private:

    //! The ID assigned
    ID _id;

    //! The name
    std::string _name;

    //! The type
    Type _type;

    //! The location of the quantity inside an element
    Location _location;

    //! The units
    std::string _units;

};


inline
SolutionDescriptor::SolutionDescriptor(const std::string& name, ID id,
    Type type, Location location, const std::string& units) :
    _name(name),
    _id(id),
    _type(type),
    _location(location),
    _units(units)
{
}


#endif /* _SOLUTIONDESCRIPTOR_H_ */
