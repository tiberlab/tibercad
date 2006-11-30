#ifndef ISE_EDGE_H_
#define ISE_EDGE_H_
#include "ISE_Vertex.h"
#include <iostream>

using  namespace  std ;

//! Edge Class. 
/*!
  Contains two Vertex pointers.
*/
class ISE_Edge
{
 public:
  //!  Constructor 
  /*! 
    Assigns two Vertex Pointers.
  */
  ISE_Edge(ISE_Vertex*  id_1,  ISE_Vertex*   id_2) ;


  //!  Virtual Destructor.
  /*!
    Dummy.
  */
  virtual ~ISE_Edge();
	
  //! Get node belonging to an  edge
  /*! Returns  a  pointer to  \c ISE_Vertex  belonging to an Edge. 
   *  v=1 -> first edge vertex; v=2 -> second edge vertex;
   */ 
  ISE_Vertex*  get_vertex(unsigned int  v); 
	
 



 private:

  /*!
    First node pointer.
  */
  ISE_Vertex*  node_1;

  /*!
    Second node pointer. 
  */
  ISE_Vertex*  node_2;

};




inline
ISE_Edge::ISE_Edge( ISE_Vertex*  id_1, ISE_Vertex*    id_2)
{
  node_1 = id_1;
  node_2 = id_2;
};

inline ISE_Vertex* 
ISE_Edge::get_vertex(unsigned int  v) 
{
  if  (v==1)
  {return node_1;}
  else if  (v==2)
  {return node_2;}
  else { cerr <<  " Error, no valid vertex in 'get_vertex' method selected "<< endl; }
	
	
};




#endif /*ISE_EDGE_H_*/
