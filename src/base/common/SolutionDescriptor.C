// $Id$

#include "SolutionDescriptor.h"

#define CASE(key, str) case SolutionDescriptor::key: os << #str; break;
#define DEFAULT default: os << "unknown"; break;


std::ostream& operator<<(std::ostream& os, SolutionDescriptor::Type type)
{
  switch (type)
  {
    CASE(REAL, real)
    CASE(COMPLEX, complex)
    CASE(VECTOR, vector)
    CASE(TENSOR, tensor)
    CASE(NTUPLE, n-tuple)
    DEFAULT
  }
  return os;
}


std::ostream& operator<<(std::ostream& os, SolutionDescriptor::Location location)
{
  switch (location)
  {
    CASE(NODES, nodal)
    CASE(CELL, cell)
    CASE(ATOM, atom)
    CASE(GLOBAL, none)
    DEFAULT
  }
  return os;
}
