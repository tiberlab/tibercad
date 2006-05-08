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
  //ISE_Element_2D(vector<ISE_Edge*> edge_id);
  //ISE_Element_2D(vector<ISE_Edge*> edge_ids);
	
  ISE_Element_2D(vector<ISE_Edge*> edge_ids, vector<bool> neg_edges  );
	
	
  virtual ~ISE_Element_2D();
  //	unsigned int get_type();
	
  //	void set_type(unsigned int element_type);
	
  vector<unsigned int> get_nodes();
	
  //	void set_nodes();
	
  //	void set_physical_region(unsigned int phys_reg);
	
 private:


  //vector<ISE_Edge> element_edges;
  vector<ISE_Edge*> element_edges;
	
  // // vector <ISE_Vertex*> element_nodes;  // in  base  class ;
	
	  
  vector<bool> negative_edges;
	
	
  void set_element_nodes();
  //	void set_element_nodes2();
	
	
  //	unsigned int  elem_type;
	
  //vector<unsigned int> element_nodes; //  in  base  class !!!
	
  //unsigned int  physical_region;
	
  void  unique_nodes_point();
	
	
  //!  Check if orientation  of  nodes  is  positive,  otherwise swap  nodes  of  element.
  /*!
   * 
   * 
   */	
  void  check_orientation_2D();
	
  void  set_element_nodes_id();

};


//inline
//ISE_Element_2D::ISE_Element_2D(vector<ISE_Edge*> edge_id ):ISE_Element() 
//{
//	for (unsigned int i=0 ;  i< edge_id.size() ; i++)
//	{
//		element_edges.push_back(*edge_id[i]);
//	}
//}

//inline
//ISE_Element_2D::ISE_Element_2D(vector<ISE_Edge*> edge_ids ):ISE_Element() 
//{
////	for (unsigned int i=0 ;  i< edge_id.size() ; i++)
//	//{
//		element_edges = edge_ids;
//		
//		element_nodes.clear();
//		set_element_nodes();
//		
//	//}
//}

inline
ISE_Element_2D::ISE_Element_2D(vector<ISE_Edge*> edge_ids, vector<bool> neg_edges  ):ISE_Element() 
{
  //	for (unsigned int i=0 ;  i< edge_id.size() ; i++)
  //{
  element_edges = edge_ids;
  negative_edges = neg_edges;
  element_nodes.clear();
  set_element_nodes();
		
  //set_element_nodes_id();
		
  // change  orientation of  element nodes if  not  counterclockwise (det > 0)
  check_orientation_2D();
		
  // writes element_nodes_id (after possible  change  of  orientation)
  set_element_nodes_id();
		
		
  //}
}




inline  vector<unsigned int> 
ISE_Element_2D::get_nodes()
{
  //  beware :  repetitions in  the  list  :  use  unique after  !
  //return  element_nodes;
  return element_nodes_id;
}


//inline unsigned int
//ISE_Element_2D::get_type()
//{
//	return elem_type;
//}

//inline void
//ISE_Element_2D::set_physical_region(unsigned int phys_reg)
//{
//	physical_region = phys_reg;
//	
//}

//inline void 
//ISE_Element_2D::set_type(unsigned int element_type)
//{
//	elem_type	= element_type;
//}



#endif /*ISE_ELEMENT_2D_H_*/
