// $Id$

#include "SpaceTransformation.h"

#include "libmesh/point.h"
#include "libmesh/tensor_value.h"


void
SpaceTransformation::rotate(const libMesh::RealVectorValue& axis,
                            const double angle,
                            libMesh::Point& point)
{
  libMesh::RealVectorValue u(axis.unit());
  double ux = u(0);
  double uy = u(1);
  double uz = u(2);

  double cost = std::cos(angle);
  double sint = std::sin(angle);

  libMesh::RealTensorValue R;

  R(0,0) = cost + ux*ux*(1-cost);
  R(0,1) = ux*uy*(1-cost) - uz*sint;
  R(0,2) = ux*uz*(1-cost) + uy*sint;
  R(1,0) = ux*uy*(1-cost) + uz*sint;
  R(1,1) = cost + uy*uy*(1-cost);
  R(1,2) = uy*uz*(1-cost) - ux*sint;
  R(2,0) = ux*uz*(1-cost) - uy*sint;
  R(2,1) = uz*uy*(1-cost) + ux*sint;
  R(2,2) = cost + uz*uz*(1-cost);

  point = R*point;
}
