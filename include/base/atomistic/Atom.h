#ifndef _ATOM_H_
#define _ATOM_H_

#include "tensor.h"
#include "TypeDefs.h"
 

//! Contains Atom definition
/*!
 * Atom is defined mainly by atomic specie and 
 * spatial vector giving the position (from library 
 * tensor.h)
 */
 class Atom
 {
	 public:
	
	 //! Atom constructor
	 Atom();
	 
	 //! Constructor with specie and position initializations
	 Atom(std::string& init_specie, Tensor1& init_position);
	 
	 //! Atom destructor
	 ~Atom();
	 	  
	 //! Atomic specie (short name)
        std::string specie;
	 
	 //! Atom position
	 Tensor1 position;
	 
	 //! An integer which says if an atom belongs to device (0)
	 //! or to contact (number of contact). Useful in electronic transport
	 unsigned int contact;

   //! ID of region containing the atom
   ID id;

   //! A general purpose integer flag (for example used in passivation)
   unsigned int flag;
	 
 };


#endif // _ATOM_H_

