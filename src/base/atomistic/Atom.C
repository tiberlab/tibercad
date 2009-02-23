#include "Atom.h"


Atom::Atom()
:_specie("none"),
belong_to_structure(false)
{
  _region_id = 0;
  _flag = 0;
  _atom_id = 0;
  _position(1) = 0.0; _position(2) = 0.0; _position(3) = 0.0;
  _conv_address[0] = 0; _conv_address[1] = 0; _conv_address[3] = 0;
  _contact = 0;
}

Atom::Atom(std::string& specie, Tensor1& position)
:belong_to_structure(false)
{
  _position = position;
  _specie = specie;
  _region_id = 0;
  _flag = 0;
  _atom_id = 0;
  _conv_address[0] = 0; _conv_address[1] = 0; _conv_address[3] = 0;
  _contact = 0;
}

Atom::Atom(std::string& specie, Tensor1& position, int (&conv_address)[3], ID atom_id, ID region_id, ID contact, unsigned int flag)
:belong_to_structure(false)
{
  _position = position;
  _specie = specie;
  _region_id = region_id;
  _flag = flag;
  _atom_id = atom_id;
  _conv_address[0] = conv_address[0]; _conv_address[1] = conv_address[1]; _conv_address[2] = conv_address[2];
  _contact = contact;
}

Atom::~Atom(){
}




