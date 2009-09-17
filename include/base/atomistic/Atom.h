#ifndef _ATOM_H_
#define _ATOM_H_

#include "tensor.h"
#include "TypeDefs.h"

class elem;


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

  //!Complete constructor: specifies all atom characteristics
  Atom(std::string& specie, Tensor1& position, int (&conv_address)[3], ID atom_id, ID region_id, unsigned int flag, ID contact);

  //! Atom destructor
  ~Atom();

  //! Set atom specie name
  void set_specie(const std::string& sp);

  //! Get atom specie name
  const std::string& get_specie() const;

  //!Set atom position (1X3 Tensor is used)
  void set_position(const Tensor1 pos);

  //! Get position coordinate i (x=1, y=2, z=3)
  const double get_position(int i) const;

  //! Get the whole position (1X3 Tensor)
  const Tensor1& get_position() const;

  //! If atom belongs to contact, return the id of contact (0 otherwise)
  const ID get_contact() const;

  //! Get the ID of the region the atom belongs to
  //! (Note: little errors may occur using basis or conventional cell
  //! preservation in Atomistic Generator)
  const int get_region_ID() const;

  //! Get atom identifier
  const ID get_atom_ID(void) const;

  //! Set atom identifier
  void set_atom_ID(const ID my_id);

  //! Set a general purpose integer flag, used internally. 0 is default safe value
  void set_flag(const unsigned int fg);

  //! Get the general purpose flag
  const unsigned int get_flag() const;

  //! True if atom belong to structure. Useful during structure construction
  bool belong_to_structure;

  //!Set element
  void set_elem(Elem* const el);

  //!Get element
  const Elem* get_elem() const;

private:

  //!Element atom belongs to
  Elem* _el;

  //! Atomic specie (short name)
  std::string _specie;

  //! Atom position
  Tensor1 _position;

  //! An integer which says if an atom belongs to device (0)
  //! or to contact (number of contact). Useful in electronic transport
  ID _contact;

  //! A general purpose integer flag (for example used in passivation)
  unsigned int _flag;

  //!An ID identifying univocally the atom
  ID _atom_id;

};


//----------------------------------------------------
// Inline member functions
//----------------------------------------------------

inline
void Atom::set_specie(const std::string& sp)
{
  _specie=sp;
}


inline
const std::string& Atom::get_specie(void) const
{
  return _specie;
}


inline
void Atom::set_position(const Tensor1 pos)
{
  _position = pos;
}


inline
const Tensor1& Atom::get_position(void) const
{
  return _position;
}


inline
const double Atom::get_position(int i) const
{
  return _position(i);
}


inline
const ID Atom::get_contact(void) const
{
  return _contact;
}


inline
const ID Atom::get_atom_ID(void) const
{
  return _atom_id;
}


inline
void Atom::set_atom_ID(const ID my_id)
{
  _atom_id=my_id;
}


inline
void Atom::set_flag(const unsigned int fg)
{
  _flag = fg;
}


inline
const unsigned int Atom::get_flag() const
{
  return _flag;
}


inline
const Elem* Atom::get_elem() const
{
  return _el;
}

#endif // _ATOM_H_

