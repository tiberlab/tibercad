#ifndef ISE_ELEMENT_H_
#define ISE_ELEMENT_H_
#include <vector>
#include "ISE_Vertex.h"

using  namespace std;

class ISE_Element
{
 public:
  ISE_Element(void);
  virtual ~ISE_Element();
	
  //virtual unsigned int get_type() = 0;
  unsigned int get_type();
	
  void set_type(unsigned int element_type);
  void set_dimension(unsigned int element_dim);
  unsigned int get_dimension();
  void set_physical_region(unsigned int phys_reg);
  unsigned int get_physical_region();

 protected:
  unsigned int  elem_type;
  unsigned int  elem_dimension;
  unsigned int  physical_region;		
		
  //	vector <unsigned int> element_nodes;
		
  vector <ISE_Vertex*> element_nodes;
  vector <unsigned int> element_nodes_id;
		
		
 private:	
  //	unsigned int  elem_type;
	
	
};

inline unsigned int 
ISE_Element::get_type() 
{
  return elem_type;
}

inline void 
ISE_Element::set_type(unsigned int element_type)
{
  elem_type	= element_type;
}

inline  unsigned int 
ISE_Element::get_dimension()
{
  return elem_dimension  ; 
}


inline void 
ISE_Element::set_dimension(unsigned int element_dim)
{
  elem_dimension = element_dim ; 
}

inline void
ISE_Element::set_physical_region(unsigned int phys_reg)
{
  physical_region = phys_reg;
	
}

inline unsigned int
ISE_Element::get_physical_region()
{
  return physical_region ;
	
}


#endif /*ISE_ELEMENT_H_*/
