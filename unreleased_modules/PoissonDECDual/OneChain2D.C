#include "OneChain2D.h"


OneChain2D::OneChain2D(const libMesh::Elem *elem, unsigned int side)
: _elem(elem),
  _side(side)
{
}
