// $Id$

#include "Scaling.h"

Scaling::Scaling(void)
  : _type(NONE),
    _is_inverted(false),
    _potential(1.0),
    _length(1.0),
    _mobility(1.0),
    _density(1.0)
{
}

Scaling::Scaling(const Scaling& scaling)
  : _type(scaling._type),
    _is_inverted(scaling._is_inverted),
    _potential(scaling._potential),
    _length(scaling._length),
    _mobility(scaling._mobility),
    _density(scaling._density)
{
}

void
Scaling::invert(void)
{
  _is_inverted = !_is_inverted;
  _potential = 1.0 / _potential;
  _length = 1.0 / _length;
  _mobility = 1.0 / _mobility;
  _density = 1.0 / _density;
}
