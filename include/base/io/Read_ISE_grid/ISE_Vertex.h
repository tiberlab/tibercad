#ifndef ISE_VERTEX_H_
#define ISE_VERTEX_H_

#include <vector>

using  namespace std;

class ISE_Vertex
{
 public:
  /*!
   *  Constructor
   */
 
  //ISE_Vertex(double  x ,  double  y ,  double z, unsigned int i);
	
  ISE_Vertex( vector <double> &  node_coord , unsigned int i );
  virtual ~ISE_Vertex();
	
	
  vector<double> get_coord();
  unsigned int get_node_id();

 private:

  vector<double> node_coordinates ;
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
