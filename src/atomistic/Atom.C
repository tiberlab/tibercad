/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file Atom.C
 * \brief tiberCAD API implementation.
 */

#include "tibercad/atomistic/Atom.h"
#include "elem.h"


Atom::Atom()
:_specie(),
_el(NULL),
_position()
{
  _virtual_type = ILLEGAL_VALUE;
}

Atom::Atom(std::string& specie, Tensor1& position)
:_el(NULL),
_position(position(1), position(2), position(3)),
_specie(specie)
{
  _virtual_type = ILLEGAL_VALUE;
  _label = ILLEGAL_VALUE;
}

Atom::Atom(std::string& specie, Tensor1& position, label_t label, atom_t type)
:_el(NULL),
_position(position(1), position(2), position(3)),
 _specie(specie),
 _label(label),
 _virtual_type(type)
{
}


void Atom::set_elem(const libMesh::Elem* el)
{
  _el = el;
}


int Atom::get_region_ID(void) const
{
  if (_el == NULL) return INVALID_ID;
  else return _el->subdomain_id();
}


void Atom::set_position(const Tensor1& pos)
{
  _position(0) = pos(1); _position(1) = pos(2); _position(2) = pos(3);
}


void
Atom::set_position(int i, double x)
{
  _position(i) = x;
}

Tensor1 Atom::get_ttype_position(void) const
{
  Tensor1 pos;
  pos(1) = _position(0); pos(2) = _position(1); pos(3) = _position(2);
  return pos;
}
