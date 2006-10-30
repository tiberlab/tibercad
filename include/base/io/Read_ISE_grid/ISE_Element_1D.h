#ifndef ISE_ELEMENT_1D_H_
#define ISE_ELEMENT_1D_H_

#include "ISE_Element.h"
#include "ISE_Vertex.h"

class ISE_Element_1D: public  ISE_Element
{
 public:

  ISE_Element_1D(ISE_Vertex*  vertex_0,  ISE_Vertex*  vertex_1);
	
  virtual ~ISE_Element_1D();
	
  vector<unsigned int> get_nodes_id();
	
	
 private:

  ISE_Vertex*  node_1;
  ISE_Vertex*  node_2;
	
};

inline 
ISE_Element_1D::ISE_Element_1D(ISE_Vertex*  vertex_0, ISE_Vertex*  vertex_1 ):ISE_Element()
{
  node_1 = vertex_0;
  node_2 = vertex_1;
		
}

inline vector<unsigned int>
ISE_Element_1D::get_nodes_id()
{
		
  element_nodes_id.push_back (node_1->get_node_id() );
  element_nodes_id.push_back (node_2->get_node_id() );
	
  return element_nodes_id;
}


#endif /*ISE_ELEMENT_1D_H_*/
