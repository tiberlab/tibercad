// $Id$

#ifndef _NODEOBJECT_H_
#define _NODEOBJECT_H_


#include "PhysicalObject.h"

class Material;


//! Description of a node, i.e. (n-3)-D object.
class NodeObject : public PhysicalObject
{

  public:

    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    ~NodeObject(void) {};


    //! Create an edge object
    /*!
     * \param options options for this boundary
     */
    static NodeObject* create(const ModelOptions& options) TBDLLOCAL;


  protected:

    //! Construct an edge object
    NodeObject(const ModelOptions& options) : PhysicalObject(NODE, options) {};


    //\copydoc PhysicalObject::do_init()
    //void do_init(void) {};


  private:

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------




#endif /* _NODEOBJECT_H_ */
