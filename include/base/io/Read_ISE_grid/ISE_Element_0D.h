#ifndef ISE_ELEMENT_0D_H_
#define ISE_ELEMENT_0D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Vertex.h"
#include <algorithm>


using namespace  std;

class ISE_Element_0D : public ISE_Element
{
 public:

  ISE_Element_0D(ISE_Vertex* vertex_0);
	
  virtual ~ISE_Element_0D();
	
  vector<unsigned int> get_nodes_id();
	
	
 private:
  ISE_Vertex* node;
	
};

inline
ISE_Element_0D::ISE_Element_0D(ISE_Vertex* vertex_0) :ISE_Element()
{
  node = vertex_0;
		
};

inline vector<unsigned int>
ISE_Element_0D::get_nodes_id()
{
		
  element_nodes_id.push_back (node->get_node_id() );
	
  return element_nodes_id;
};


#endif /*ISE_ELEMENT_0D_H_*/
