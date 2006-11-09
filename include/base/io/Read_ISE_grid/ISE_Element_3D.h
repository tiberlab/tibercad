#ifndef ISE_ELEMENT_3D_H_
#define ISE_ELEMENT_3D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Edge.h"
#include "ISE_Vertex.h"
#include "ISE_Face.h"
#include <algorithm>

using namespace  std;

class ISE_Element_3D : public ISE_Element
{
 public:

  ISE_Element_3D(vector<ISE_Face*> face_ids, vector<bool> neg_faces  );
	
	
  virtual ~ISE_Element_3D();
	
  vector<unsigned int> get_nodes_id();
	
	
 private:


  vector<ISE_Face*> element_faces;
	  
  vector<bool> negative_faces;
	
	
  void set_element_nodes();


  void  unique_nodes_point();

	
  void  set_element_nodes_id();

};


inline
ISE_Element_3D::ISE_Element_3D(vector<ISE_Face*> face_ids, vector<bool> neg_faces  ):ISE_Element()
{
 
  element_faces = face_ids;
  negative_faces = neg_faces;
  element_nodes.clear();
  set_element_nodes();
  unique_nodes_point();
  // writes element_nodes_id (after possible  change  of  orientation)
  set_element_nodes_id();
		
}




inline  vector<unsigned int> 
ISE_Element_3D::get_nodes_id()
{
  //  beware :  repetitions in  the  list  :  use  unique after  !
  return element_nodes_id;
}


#endif /*ISE_ELEMENT_3D_H_*/
