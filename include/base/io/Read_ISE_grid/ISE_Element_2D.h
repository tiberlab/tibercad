#ifndef ISE_ELEMENT_2D_H_
#define ISE_ELEMENT_2D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Edge.h"
#include "ISE_Vertex.h"
#include <algorithm>

using namespace  std;

class ISE_Element_2D : public  ISE_Element
{
 public:
	
  ISE_Element_2D(vector<ISE_Edge*> edge_ids, vector<bool> neg_edges);
	
  virtual ~ISE_Element_2D();

  vector<unsigned int> get_nodes_id();
	
 private:


  vector<ISE_Edge*> element_edges;
	
	  
  vector<bool> negative_edges;
	


  //!  Check if orientation  of  nodes  is  positive,  otherwise swap  nodes  of  element.	
  void set_element_nodes();

  void  unique_nodes_point();
	
  void  set_element_nodes_id();

};



inline
ISE_Element_2D::ISE_Element_2D(vector<ISE_Edge*> edge_ids, vector<bool> neg_edges  ):ISE_Element()
{
  element_edges = edge_ids;
  negative_edges = neg_edges;
  element_nodes.clear();
  set_element_nodes();
  unique_nodes_point();
		
  // writes element_nodes_id (after possible  change  of  orientation)
  set_element_nodes_id();
		
		
}




inline  vector<unsigned int> 
ISE_Element_2D::get_nodes_id()
{
  //  beware :  repetitions in  the  list  :  to use after unique_nodes_point() method!
  return element_nodes_id;
}




#endif /*ISE_ELEMENT_2D_H_*/
