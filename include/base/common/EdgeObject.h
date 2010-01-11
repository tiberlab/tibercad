// $Id$

#ifndef _EDGEOBJECT_H_
#define _EDGEOBJECT_H_


#include "PhysicalObject.h"

class Material;


//! Description of an edge, i.e. (n-2)-D object.
class EdgeObject : public PhysicalObject
{

  public:

    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    ~EdgeObject(void) {};


    //! Create an edge object
    /*!
     * \param options options for this boundary
     */
    static EdgeObject* create(const ModelOptions& options);


  protected:

    //! Construct an edge object
    EdgeObject(void) : PhysicalObject(EDGE) {};


    //! \copydoc PhysicalObject::do_init()
    void do_init(void);


  private:

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
EdgeObject*
EdgeObject::create(const ModelOptions& options)
{
  return new EdgeObject();
}



#endif /* _EDGEOBJECT_H_ */
