#ifndef ISE_VERTEX_H_
#define ISE_VERTEX_H_

#include "tiber_dll.h"

#include <vector>

using  namespace std;

//! Vertex Class.
/*! 
  Contains node coordinates and id.
*/
class TBDLLOCAL ISE_Vertex
{
 public:
  //! Constructor.
  /*!
    Assigns node coordinates and id.
  */
  ISE_Vertex( vector <double> &  node_coord , unsigned int i );

  //!Virtual Destructor.
  /*! 
    Dummy.
  */
  virtual ~ISE_Vertex();
	
  /*!
    Returns coordinates vector.
  */	
  vector<double> get_coord();

  /*!
    Returns node id.
  */
  unsigned int get_node_id(); 

 private:

  /*!
    Vertex coordinates.
  */
  vector<double> node_coordinates ;

  /*!
    Node id.
  */
  unsigned int node_ID;	
	
};

inline 
ISE_Vertex::ISE_Vertex(vector<double>&  node_coord , unsigned int i) 
{
  node_coordinates = node_coord;
  node_ID = i;
	
}


inline vector<double> 
ISE_Vertex::get_coord()
{
  return node_coordinates;
}

inline unsigned int 
ISE_Vertex::get_node_id()
{
  return node_ID;
}


#endif /*ISE_VERTEX_H_*/
