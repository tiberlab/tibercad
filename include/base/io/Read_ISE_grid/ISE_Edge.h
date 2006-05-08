#ifndef ISE_EDGE_H_
#define ISE_EDGE_H_
#include "ISE_Vertex.h"
#include <iostream>

//class ISE_Vertex;
using  namespace  std ;

class ISE_Edge
{
 public:

  // constructor
  ISE_Edge(ISE_Vertex*  id_1,  ISE_Vertex*   id_2) ;
	
  //! Get node belonging to an  edge
  /*!
   * Returns  a  pointer to  \c ISE_Vertex 
   * (v=1 -> first node of  edge,  v=2 ->  second node of  edge
   */
  ISE_Vertex*  get_vertex(unsigned int  v);
	
  virtual ~ISE_Edge();
	
 private:
  ISE_Vertex*   node_1;
  ISE_Vertex*  node_2;

};


inline
ISE_Edge::ISE_Edge( ISE_Vertex*  id_1, ISE_Vertex*    id_2)
{
  node_1 = id_1;
  node_2 = id_2;
}

inline ISE_Vertex* 
ISE_Edge::get_vertex(unsigned int  v)
{
  if  (v==1)
  {return node_1;}
  else if  (v==2)
  { return node_2;}
  else { cerr <<  " Error "<< endl; }
	
	
}



#endif /*ISE_EDGE_H_*/
