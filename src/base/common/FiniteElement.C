// $Id$

#include "FiniteElement.h"
#include "Constants.h"



#include "elem.h"
#include "quadrature.h"
#include "fe_map.h"
//#include "fe_abstract.h"
//#include "fe_base.h"
#include "libMeshDefs.h"

using namespace std;

//USELIBMESHTYPE(FE);


template <unsigned int Dim, libMesh::FEFamily T>
void
FiniteElement<Dim, T>::reinit(const libMesh::Elem* elem, const unsigned int side,
    const double tolerance, const std::vector<libMesh::Point>* const,
    const std::vector<double>* const)
{
  // we need this, otherwise below it will complain in some cases when accessing
  this->_fe_map->get_JxW();

  libMesh::FE<Dim, T>::reinit(elem, side, tolerance);

  unsigned int dim = libMesh::FE<Dim, T>::qrule->get_dim();

  double x0 = _length_scaling / _mesh_units;
  double x0_inv = _mesh_units / _length_scaling;

  double J = x0_inv;
  switch (dim)
  {
    case 3:
      J *= x0_inv;
    case 2:
      J *= x0_inv;
      break;
    case 0: // in 1D simulations!
      J = 1.0;
  }


  // TODO : check the _fe_map issue!
  unsigned int n_points = this->_fe_map->get_JxW().size();
  switch (_symmetry)
  {
    case TiberCad::CYLINDRICAL:
      for (unsigned int i = 0; i < n_points; i++){
        double r = abs(this->_fe_map->get_xyz()[i](0));
        this->_fe_map->get_JxW()[i] *= 2 * M_PI * r * x0_inv * J;
      }
      break;
    default:
      for (unsigned int i = 0; i < n_points; i++)
        this->_fe_map->get_JxW()[i] *= J;
      break;
  }

  // is already done because FE<Dim, T>::reinit(elem, points) gets called
  // automatically
  //if (FE<Dim, T>::calculate_dphi)
  //  for (unsigned int i = 0; i < n_points; i++)
  //    for (unsigned int j = 0; j < elem->n_nodes(); j++)
  //      FE<Dim, T>::dphi[j][i] *= x0;

}


template <unsigned int Dim, libMesh::FEFamily T>
void
FiniteElement<Dim, T>::reinit(const libMesh::Elem* elem, const vector<libMesh::Point>* const points,
    const std::vector<double>* const)
{
  // we need this, otherwise below it will complain in some cases when accessing
  this->_fe_map->get_JxW();

  libMesh::FE<Dim, T>::reinit(elem, points);

  unsigned int dim = libMesh::FE<Dim, T>::dim;

  double x0 = _length_scaling / _mesh_units;
  double x0_inv = _mesh_units / _length_scaling;

  double J = x0_inv;
  switch (dim)
  {
    case 3:
      J *= x0_inv;
    case 2:
      J *= x0_inv;
      break;
  }


  unsigned int n_points = this->_fe_map->get_JxW().size();
  switch (_symmetry)
  {
    case TiberCad::CYLINDRICAL:
      for (unsigned int i = 0; i < n_points; i++){
        double r = abs(this->_fe_map->get_xyz()[i](0));
        this->_fe_map->get_JxW()[i] *= 2 * M_PI * r * x0_inv * J;
        //libMesh::FEBase::JxW[i] *= 2 * M_PI * abs(libMesh::FEBase::xyz[i](0)) * x0_inv * J;
      }
      break;
    default:
      for (unsigned int i = 0; i < n_points; i++)
      {
        this->_fe_map->get_JxW()[i] *= J;
      }
      break;
  }

  if (libMesh::FE<Dim, T>::calculate_dphi)
  {
    unsigned int n_dofs = libMesh::FE<Dim, T>::dphi.size();
    for (unsigned int i = 0; i < n_points; i++)
      for (unsigned int j = 0; j < n_dofs; j++)
        libMesh::FE<Dim, T>::dphi[j][i] *= x0;
  }
}





// explicit instantiations

template class FiniteElement<1, libMeshEnums::LAGRANGE>;
template class FiniteElement<2, libMeshEnums::LAGRANGE>;
template class FiniteElement<3, libMeshEnums::LAGRANGE>;

template class FiniteElement<1, libMeshEnums::MONOMIAL>;
template class FiniteElement<2, libMeshEnums::MONOMIAL>;
template class FiniteElement<3, libMeshEnums::MONOMIAL>;
