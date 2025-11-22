// $Id$

#include "tibercad/math/Scaling.h"

Scaling::Scaling(void)
  : _type(NONE),
    _potential(1.0),
    _length(1.0),
    _mesh_units(1.0),
    _mobility(1.0),
    _density(1.0)
{
}

Scaling::Scaling(const Scaling& scaling)
  : _type(scaling._type),
    _potential(scaling._potential),
    _length(scaling._length),
    _mesh_units(scaling._mesh_units),
    _mobility(scaling._mobility),
    _density(scaling._density)
{
}

