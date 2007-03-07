#ifndef ISE_ELEMENT_2D_H_
#define ISE_ELEMENT_2D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Edge.h"
#include "ISE_Vertex.h"
#include <algorithm>

using namespace  std;

//!2D Element Class.
/*!
  Contains Segment data and orientations.
*/
class ISE_Element_2D : public  ISE_Element
{
 public:

  //!  Constructor
  /*!
    Assigns Edge pointers and defines orientation vector.
  */	
  ISE_Element_2D(vector<ISE_Edge*> edge_ids, vector<bool> neg_edges);

  //!  Virtual Destructor.
  /*!
    Dummy.
  */	
  virtual ~ISE_Element_2D();

  /*!
   *  Returns nodes id vector.
   *  Beware: repetitions in the list: to use after unique_nodes_point() method!
   */
  vector<unsigned int> get_nodes_id();
	
 private:

  /*!
    Edge Pointers vector.
  */
  vector<ISE_Edge*> element_edges;
	
  /*!
    Orientation values vector. A 'false' value determines opposite orientation.
  */	  
  vector<bool> negative_edges;
	

  /*!
    Check if orientation  of  nodes  is  positive,  otherwise swap  nodes  of  element.
  */	
  void set_element_nodes();


  /*!
    Unique Node Verification. Deletes repetitions in the list.
  */
  void  unique_nodes_point();


  /*!
    Writes element_nodes_id vector. Used after possible change of orientation.
  */	
  void  set_element_nodes_id();

  /*!
    Check orientation of the element and change order of  nodes if  necessary.
  */
  void  check_orientation_2D();

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
