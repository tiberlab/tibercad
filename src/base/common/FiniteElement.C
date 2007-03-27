// $Id$

#include "FiniteElement.h"
#include "quadrature.h"

#include "elem.h"

using namespace std;



template <unsigned int Dim, FEFamily T>
void
FiniteElement<Dim, T>::reinit(const Elem* elem, const unsigned int side)
{
  FE<Dim, T>::reinit(elem, side);

  unsigned int dim = FE<Dim, T>::qrule->get_dim();
  unsigned int n_points = FE<Dim, T>::JxW.size();

  double x0_inv = _mesh_units / _length_scaling;

  double J = x0_inv;
  switch (dim)
  {
    case 3:
      J *= x0_inv;
    case 2:
      J *= x0_inv;
  }

  switch (_symmetry)
  {
    case TiberCad::CYLINDRICAL:
      for (unsigned int i = 0; i < n_points; i++)
        FE<Dim, T>::JxW[i] *= 2 * M_PI * FE<Dim, T>::xyz[i](0) * x0_inv * J;
      break;
    default:
      for (unsigned int i = 0; i < n_points; i++)
        FE<Dim, T>::JxW[i] *= J;
      break;
  }

}


template <unsigned int Dim, FEFamily T>
void
FiniteElement<Dim, T>::reinit(const Elem* elem, const vector<Point>* points)
{
  FE<Dim, T>::reinit(elem, points);

  unsigned int dim = FE<Dim, T>::dim;
  unsigned int n_points = FE<Dim, T>::JxW.size();

  double x0 = _length_scaling / _mesh_units;
  double x0_inv = _mesh_units / _length_scaling;

  double J = x0_inv;
  switch (dim)
  {
    case 3:
      J *= x0_inv;
    case 2:
      J *= x0_inv;
  }

  switch (_symmetry)
  {
    case TiberCad::CYLINDRICAL:
      for (unsigned int i = 0; i < n_points; i++)
        FE<Dim, T>::JxW[i] *= 2 * M_PI * FE<Dim, T>::xyz[i](0) * x0_inv * J;
      break;
    default:
      for (unsigned int i = 0; i < n_points; i++)
        FE<Dim, T>::JxW[i] *= J;
      break;
  }
  
  //if (FE<Dim, T>::calculate_dphi)
    for (unsigned int i = 0; i < n_points; i++)
      for (unsigned int j = 0; j < elem->n_nodes(); j++)
        FE<Dim, T>::dphi[j][i] *= x0;
}





// explicit instantiations

template class FiniteElement<1, libMeshEnums::LAGRANGE>;
template class FiniteElement<2, libMeshEnums::LAGRANGE>;
template class FiniteElement<3, libMeshEnums::LAGRANGE>;
