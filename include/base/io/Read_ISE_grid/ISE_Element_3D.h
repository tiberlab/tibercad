#ifndef ISE_ELEMENT_3D_H_
#define ISE_ELEMENT_3D_H_
#include <vector>
#include "ISE_Element.h"
#include "ISE_Edge.h"
#include "ISE_Vertex.h"
#include "ISE_Face.h"
#include <algorithm>

using namespace  std;

//! 3D Element Class.
/*!
  Contains faces data and orientations.
*/
class TBDLLOCAL ISE_Element_3D : public ISE_Element
{
 public:
  //!  Constructor
  /*!
    Assigns face pointers and orientations vector.
  */
  //  ISE_Element_3D(vector<ISE_Face*> face_ids, vector<bool> neg_faces  );
  ISE_Element_3D(vector<ISE_Face*> face_ids, vector<bool> neg_faces,unsigned int element_type);

  //! Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Element_3D();


 private:

  /*!
    Face pointers vector.
  */
  vector<ISE_Face*> element_faces;


  /*!
    Orientation vector. If negative, face must be inverted.
  */
  vector<bool> negative_faces;

  /*!
    Sets Element nodes.
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
  void  check_orientation_3D();

  // virtual void set_type(unsigned int element_type);

};


inline
ISE_Element_3D::ISE_Element_3D(vector<ISE_Face*> face_ids, vector<bool> neg_faces,
                               unsigned int element_type ):ISE_Element()
{

  element_faces = face_ids;
  negative_faces = neg_faces;
  set_type(element_type);
  element_nodes.clear();
  set_element_nodes();
  unique_nodes_point();
  // writes element_nodes_id (after possible  change  of  orientation)
  set_element_nodes_id();

}





/* inline void  */
/* ISE_Element_3D::set_type(unsigned int element_type) */
/* { */
/*   elem_type = element_type; */
/*   // check_orientation_3D(); */

 /* } */

#endif /*ISE_ELEMENT_3D_H_*/
