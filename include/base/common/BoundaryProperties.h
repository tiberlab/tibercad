// $Id$

#ifndef _BOUNDARYPROPERTIES_H_
#define _BOUNDARYPROPERTIES_H_


#include "ModelOptions.h"
#include "TypeDefs.h"


//! Boundary properties for a certain type of simulation.
/*!
 *  This is the base class for boundary properties. Every solver
 *  module should implement its own class derived from this one to hold
 *  the properties needed for calculations.
 */
class BoundaryProperties
{

  public:

    //! The empty destructor
    /*!
     * should be implemented in the derived classes if needed
     */
    virtual ~BoundaryProperties(void) {};

    //! Initialize this boundary
    void init(void);


    //! Get ID of the boundary
    /*!
     * The ID of a boundary should be a unique number
     */
    //ID get_id(void) const;


  protected:
    
    //! The empty constructor
    /*!
     * \c BoundaryType should not be instantiated directly
     */
    BoundaryProperties(void);
 
    //! Set options for this model
    /*!
     * The options are stored internally and are accessible through
     * special methods.
     * Options have to be specified at creation time.
     */
    void set_options(const ModelOptions& options);

    //! Get the options for this contact model
    const ModelOptions& get_options(void) const;
   
    //! Initialize the model
    /*!
     * This method should set all model options and call
     * \c init() of any associated model
     *
     * This method should be implemented in derived classes.
     */
    virtual void do_init(void) {};


  private:

    //! Disable copy constructor
    BoundaryProperties(const BoundaryProperties&);

    //! Disable assignement operator
    BoundaryProperties& operator=(const BoundaryProperties&); 

    //! The identifier of a boundary object
    /*!
     * The identifier is assumed to be uniqe for every instance of
     * \c Boundary
     */
    //ID _id;

    //! The options for this model as read from the input file
    ModelOptions _options;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
BoundaryProperties::BoundaryProperties(void)
{
}


inline
void
BoundaryProperties::init(void)
{
  do_init();
}

inline
void
BoundaryProperties::set_options(const ModelOptions& options)
{
  _options += options;
}


inline
const ModelOptions&
BoundaryProperties::get_options(void) const
{
  return _options;
}


#endif // _BOUNDARYPROPERTIES_H_
