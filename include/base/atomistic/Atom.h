// $Id$

#ifndef _ATOM_H_
#define _ATOM_H_

#include "tensor.h"
#include "TypeDefs.h"
#include "point.h"
#include "dof_object.h"
#include "Specie.h"


#include "elem.h"

#define ILLEGAL_VALUE 255 



//! Contains Atom definition
/*!
 * Atom is defined mainly by atomic specie and
 * spatial vector giving the position (from library
 * tensor.h)
 */
class Atom
{
public:

  typedef unsigned char label_t;
  typedef unsigned char atom_t;

  //! Atom constructor
  Atom();

  //! Constructor with specie and position initializations
  Atom(std::string& init_specie, Tensor1& init_position);

  //!Complete constructor: specifies all atom characteristics
  Atom(std::string& specie, Tensor1& position, label_t label, atom_t type);

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

  //! Set position
  void set_position(const libMesh::Point& p);

  /*! \brief Get absolute position coordinate i (x=0, y=1, z=2) 
   * atomic coordinates are stored in Angstrom
   */
  double get_position(int i) const;

  //! Get the whole position (1X3 Tensor)
  Tensor1 get_ttype_position() const;

  //! Get the position as a Point object
  const libMesh::Point& get_position() const;

  /*! \brief Get the ID of the region the atom belongs to
   * (Note: it will get the region from associated element)
   */
  int get_region_ID() const;

   /*! \brief Set label
    *   
    *  Label is used to mark the atom number within the primitive cell
    *  In binary materials this is equivalent to cation/anion species
    *  Database follows the convention that specie_1 is cation
    *  Valid numbers are in the range 1 to 255 (unsigned char)
    *  0 is reserved for passivation atoms.
    */
  void set_label(unsigned int fg);

  //! Get atomic label
  label_t get_label(void) const;

  //! checks whether an atom is cation, which by CONVENTION is label ==1 
  bool is_cation(void) const;
  
  //! Get the general purpose flag
  void set_type(atom_t fg);

  //! Get the general purpose flag
  atom_t get_type(void) const;


  //!Set element
  void set_elem(const libMesh::Elem* el);

  //!Get element
  const libMesh::Elem* get_elem() const;

private:

  //!Element atom belongs to
  const libMesh::Elem* _el;

  //! Atomic specie (short name)
  Specie _specie;

  //! Atom position
  libMesh::Point _position;

  //! Atom number in primitive cell (e.g. anion/cation)
  label_t _label;

  //! Used for special purposes
  atom_t _virtual_type;

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
void Atom::set_position(const libMesh::Point& p)
{
  _position = p;
}

inline 
const libMesh::Point& Atom::get_position() const
{
  return _position;
}  


inline
void Atom::set_label(unsigned int fg)
{
  _label = static_cast<Atom::label_t>(fg);
}

inline
Atom::label_t Atom::get_label(void) const
{
  return _label;
}

inline
bool Atom::is_cation(void) const
{
  return ((get_label() == 1) || (get_label() == 3));
}

inline
void Atom::set_type(Atom::atom_t fg)
{
  _virtual_type = fg;
}
/*
inline
void Atom::set_type(unsigned int fg)
{
  _virtual_type = static_cast<Atom::atom_t>(fg);
}
*/

inline
Atom::atom_t Atom::get_type() const
{
  return _virtual_type;
}




inline
const libMesh::Elem* Atom::get_elem() const
{
  return _el;
}

#endif // _ATOM_H_

