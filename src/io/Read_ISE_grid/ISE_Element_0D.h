#ifndef ISE_ELEMENT_0D_H_
#define ISE_ELEMENT_0D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Vertex.h"
#include <algorithm>


using namespace  std;


//! 0D Class Element.
/*!
  Point Definition.
*/
class TBDLLOCAL ISE_Element_0D : public ISE_Element
{
 public:
  //! Constructor.
  /*!
    Assigns a Vertex pointer to the element.
  */
  ISE_Element_0D(ISE_Vertex* vertex_0);

  //! Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Element_0D();



 private:

  /*!
    Vertex Pointer.
  */
  ISE_Vertex* node;

};

inline
ISE_Element_0D::ISE_Element_0D(ISE_Vertex* vertex_0) :ISE_Element()
{
  node = vertex_0;
  element_nodes_id.resize(1);
  element_nodes_id[0] = node->get_node_id();

}



#endif /*ISE_ELEMENT_0D_H_*/
