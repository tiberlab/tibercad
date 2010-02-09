#include "Atom.h"
#include "elem.h"
#include"elem.h"


Atom::Atom()
:_specie("none"),
belong_to_structure(false),
_el(NULL),
_position()
{
  _flag = 0;
}

Atom::Atom(std::string& specie, Tensor1& position)
:belong_to_structure(false),
_el(NULL),
_position(position(1), position(2), position(3)),
_specie(specie)
{
  _flag = 0;
}

Atom::Atom(std::string& specie, Tensor1& position, unsigned int flag)
:belong_to_structure(false),
_el(NULL),
_position(position(1), position(2), position(3)),
_specie(specie)
{
  _flag = flag;
}

Atom::~Atom()
{
}


void Atom::set_elem(Elem* el)
{
  _el = el;
}


const int Atom::get_region_ID(void) const
{
  if (_el == NULL) return INVALID_ID;
  else return _el->subdomain_id();
}


void Atom::set_position(const Tensor1 pos)
{
  _position(0) = pos(1); _position(1) = pos(2); _position(2) = pos(3);
}


Tensor1 Atom::get_position(void) const
{
  Tensor1 pos;
  pos(1) = _position(0); pos(2) = _position(1); pos(3) = _position(2);
  return pos;
}
