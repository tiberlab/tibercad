// $Id$

#ifndef _NODEOBJECT_H_
#define _NODEOBJECT_H_


#include "PhysicalObject.h"

class Material;


//! Description of an edge, i.e. (n-2)-D object.
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
    static NodeObject* create(const ModelOptions& options);


  protected:

    //! Construct an edge object
    NodeObject(const ModelOptions& options) : PhysicalObject(NODE, options) {};


    //! \copydoc PhysicalObject::do_init()
    void do_init(void);


  private:

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
NodeObject*
NodeObject::create(const ModelOptions& options)
{
  return new NodeObject(options);
}



#endif /* _NODEOBJECT_H_ */
