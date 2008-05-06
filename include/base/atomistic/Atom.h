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
	 	
	 // Meber functions to obtain private data
	 void set_specie(const std::string& sp) {specie=sp;}
	 const std::string& get_specie() const {return specie;}

	 void set_position(Tensor1 pos){position = pos;}
	 double get_position(int i) const {return position(i);}
	 Tensor1 get_position() const {return position;}

         int get_contact() const {return contact;}
	 int get_ID() const {return id;}
	 void set_ID(int my_id) {id=my_id;}

	 void set_flag(unsigned int fg){flag = fg;}
	 unsigned int get_flag() const {return flag;}



     private:

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

