#ifndef ISE_ELEMENT_1D_H_
#define ISE_ELEMENT_1D_H_

#include "ISE_Element.h"
#include "ISE_Vertex.h"


//! 1D Element Class.
/*!
  Segment Definition.
*/
class ISE_Element_1D: public  ISE_Element
{
 public:

  //!  Constructor
  /*!
    Assigns two vertex pointers to the element.
  */
  ISE_Element_1D(ISE_Vertex*  vertex_0,  ISE_Vertex*  vertex_1);

  //!Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Element_1D();



 private:

  /*!
    Vertex Pointer.
  */
  ISE_Vertex*  node_1;

  /*!
    Vertex Pointer.
  */
  ISE_Vertex*  node_2;

};

inline
ISE_Element_1D::ISE_Element_1D(ISE_Vertex*  vertex_0, ISE_Vertex*  vertex_1 ):ISE_Element()
{
  node_1 = vertex_0;
  node_2 = vertex_1;

  element_nodes_id.resize(2);
  element_nodes_id[0] = node_1->get_node_id();
  element_nodes_id[1] = node_2->get_node_id();

}



#endif /*ISE_ELEMENT_1D_H_*/
