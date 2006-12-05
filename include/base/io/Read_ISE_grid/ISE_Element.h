#ifndef ISE_ELEMENT_H_
#define ISE_ELEMENT_H_
#include <vector>
#include "ISE_Vertex.h"
#include "reference_counted_object.h"

using  namespace std;


//! Element Base Class.
/*!
  Contains general methods. 
*/
class ISE_Element : public ReferenceCountedObject<ISE_Element>
{
 public:

  //!  Constructor.
  /*!
    Dummy. 
  */
  ISE_Element(void);

  //!Virtual Destructor.
  /*!
    Dummy. 
  */
  virtual ~ISE_Element();

  /*!
    Returns Element type.
  */	
  unsigned int get_type();
	
  /*!
    Sets Element type.
  */
  void set_type(unsigned int element_type);

  /*!
    Sets Element dimension. 
  */
  void set_dimension(unsigned int element_dim);

  /*!
    Gets Element dimension.
  */
  unsigned int get_dimension();

  /*!
    Sets Physical region.
  */
  void set_physical_region(unsigned int phys_reg);

  /*!
    Gets Physical region.
  */
  unsigned int get_physical_region();


 protected:

  /*!
    Defines Element type.
  */
  unsigned int  elem_type;

  /*!
    Defines Element dimension.
  */
  unsigned int  elem_dimension;

  /*!
    Defines Physical region.
  */
  unsigned int  physical_region;		
		
  /*!
    Element Vertex pointers vector.
  */
  vector <ISE_Vertex*> element_nodes;

  /*!
    Element nodes id vector.
  */
  vector <unsigned int> element_nodes_id;
	
};

inline unsigned int 
ISE_Element::get_type() 
{
  return elem_type;
}

inline void 
ISE_Element::set_type(unsigned int element_type)
{
  elem_type = element_type;
}

inline  unsigned int 
ISE_Element::get_dimension()
{
  return elem_dimension; 
}


inline void 
ISE_Element::set_dimension(unsigned int element_dim)
{
  elem_dimension = element_dim; 
}

inline void
ISE_Element::set_physical_region(unsigned int phys_reg)
{
  physical_region = phys_reg;
	
}

inline unsigned int
ISE_Element::get_physical_region()
{
  return physical_region;
	
}


#endif /*ISE_ELEMENT_H_*/
