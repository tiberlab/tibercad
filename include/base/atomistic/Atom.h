// $Id$

#ifndef _ATOM_H_
#define _ATOM_H_

#include "tensor.h"
#include "TypeDefs.h"
#include "point.h"
#include "dof_object.h"
#include "Specie.h"

class Elem;


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
  Atom(std::string& specie, Tensor1& position, unsigned int flag);

  //! Atom destructor
  ~Atom();

  //! Set atom specie name
  void set_specie(const std::string& sp);

   //! Set atom specie name
  void set_specie(const Specie& sp);

  //! Get atom specie name
  const Specie& get_specie() const;

  //!Set atom position (1X3 Tensor is used)
  void set_position(const Tensor1 pos);

  //! Set position number i
  void set_position(int i, double x);

  //! Get absolute position coordinate i (x=0, y=1, z=2) 
  //! atomic coordinates are stored in Angstrom
  double get_position(int i) const;

  //! Get the whole position (1X3 Tensor)
  Tensor1 get_ttype_position() const;

  //! Get the position as a Point object
  Point get_position() const;

  //! Get the ID of the region the atom belongs to
  //! (Note: little errors may occur using basis or conventional cell
  //! preservation in Atomistic Generator)
  int get_region_ID() const;

  //! Set a general purpose integer flag, used internally. 0 is default safe value
  void set_flag(const unsigned int fg);

  //! Get the general purpose flag
  unsigned int get_flag() const;

  //! True if atom belong to structure. Useful during structure construction
  bool belong_to_structure;

  //!Set element
  void set_elem(Elem* el);

  //!Get element
  const Elem* get_elem() const;

private:

  //!Element atom belongs to
  Elem* _el;

  //! Atomic specie (short name)
  Specie _specie;

  //! Atom position
  Point _position;

  //! A general purpose integer flag (for example used in passivation)
  unsigned int _flag;

  Tensor1 _ttype_position;

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
void Atom::set_specie(const Specie& sp)
{
 _specie = sp;
}

inline
const Specie& Atom::get_specie(void) const
{
  return _specie;
}


inline
double Atom::get_position(int i) const
{
  return _position(i);
}


inline 
Point Atom::get_position() const
{
  return _position;
}  

inline
void Atom::set_flag(const unsigned int fg)
{
  _flag = fg;
}


inline
unsigned int Atom::get_flag() const
{
  return _flag;
}


inline
const Elem* Atom::get_elem() const
{
  return _el;
}

#endif // _ATOM_H_

