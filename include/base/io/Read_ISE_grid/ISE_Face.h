#ifndef ISE_FACE_H_
#define ISE_FACE_H_

#include "ISE_Vertex.h"
#include "ISE_Edge.h"
#include <iostream>


using  namespace  std ;

class ISE_Face
{
 public:

  ISE_Face(vector <ISE_Edge*> fcs_edg, vector <bool> neg_edg) ;
	
  //! Get edge belonging to a face
  //Returns  a  pointer to  \c ISE_Edge 
  //(e=0 -> first edge of face,  e=1 ->  second edge of face, ...
  //An overflow control is included (the "e" value is limited)
  ISE_Edge* get_edge (unsigned int e);

//Returns the state of the 'e' edge (true == negative edge, false == positive edge)
  bool get_edge_state (unsigned int e);

//Returns the number of edges in the face
  unsigned int get_face_edges_size();
	
  virtual ~ISE_Face();
	
 private:

  vector <ISE_Edge*> face_edges;
  vector <bool> negative_face_edges;
};


inline ISE_Face::ISE_Face(vector <ISE_Edge*> fcs_edg, vector <bool> neg_edg)
{
	
	for (unsigned int e=0; e < fcs_edg.size(); e++)
	{
		face_edges.push_back(fcs_edg[e]);
		negative_face_edges.push_back(neg_edg[e]);
	}

};



inline ISE_Edge* ISE_Face::get_edge (unsigned int e)
{

	if (e < face_edges.size())
	{
		return face_edges[e];
	}
	else
	{
		cerr << "The 'e' value for the 'ISE_Edge* ISE_Face::get_edge(unsigned int e)' method is incorrect" << endl;
		exit(1);
	}
};


inline bool ISE_Face::get_edge_state (unsigned int e)
{

	if (e < negative_face_edges.size())
	{
		return negative_face_edges[e];
	}
	else
	{
		cerr << "The 'e' value for the 'bool ISE_Face::get_edge_state(unsigned int e)' method is incorrect" << endl;
		exit(1);
	}
};



inline unsigned int ISE_Face::get_face_edges_size()
{
	return face_edges.size();
};

#endif /*ISE_FACE_H_*/
