#include "Atom.h"
#include "elem.h"
#include"elem.h"


Atom::Atom()
:_specie("none"),
belong_to_structure(false),
_el(NULL)
{
  _flag = 0;
  _atom_id = 0;
  _position(1) = 0.0; _position(2) = 0.0; _position(3) = 0.0;
  _contact = 0;
}

Atom::Atom(std::string& specie, Tensor1& position)
:belong_to_structure(false),
_el(NULL)
{
  _position = position;
  _specie = specie;
  _flag = 0;
  _atom_id = 0;
  _contact = 0;
}

Atom::Atom(std::string& specie, Tensor1& position, int (&conv_address)[3], ID atom_id, ID region_id, ID contact, unsigned int flag)
:belong_to_structure(false),
_el(NULL)
{
  _position = position;
  _specie = specie;
  _flag = flag;
  _atom_id = atom_id;
  _contact = contact;
}

Atom::~Atom()
{
}


void Atom::set_elem(Elem* const el)
{
  _el = el;
}


const int Atom::get_region_ID(void) const
{
  if (_el == NULL) return INVALID_ID;
  else return _el->subdomain_id();
}
