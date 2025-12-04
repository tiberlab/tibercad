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
 * \file Kspace.C
 * \brief tiberCAD API implementation.
 */


#include "tibercad/kintegration/Kspace.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/physics/Constants.h"
#include "tibercad/io/Messages.h"
#include "tibercad/utils/Utils.h"
#include "tibercad/math/SpaceTransformation.h"
#include "tibercad/math/TensorOperators.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/cell_hex8.h"

#include "libmesh/mesh_tetgen_interface.h"


#include <cmath>



using namespace std;
using namespace libMesh;


namespace
{
  bool compare_points(const Point& a, const Point& b)
  {
    return(a.absolute_fuzzy_equals(b));
  }
}

std::vector<std::string>
Kspace::symmetry_names = {"Gamma", "linear", "quadratic", "rectangular",
                          "hexagonal", "fcc", "bcc", "cubic",
                          "tetragonal", "orthorhombic"};

//------------------------------------------------------------------------------//
Kspace::Kspace(const ModelOptions& options, const libMesh::Parallel::Communicator& comm)
 : mod_opt(options),
   _mesh_order(libMesh::FIRST),
   kmesh(nullptr),
   k_max(1, 1, 1),
   k_space_symmetry(GAMMA)
{

  // we create a serial communicator
  comm.split(0, 0, kspace_comm); 

  transform_matrix = Tensor2(1);

  do_init();
}


// Copy constructor
Kspace::Kspace( const Kspace& kspace)
 : mod_opt(kspace.mod_opt),
   degeneracy_factor(kspace.degeneracy_factor),
   _mesh_order(kspace._mesh_order),
   k_space_dim(kspace.k_space_dim),
   k_basis_vector1(kspace.k_basis_vector1),
   k_basis_vector2(kspace.k_basis_vector2),
   k_basis_vector3(kspace.k_basis_vector3),
   b1(0),
   b2(1),
   b3(2),
   k_max(kspace.k_max),
   k_space_symmetry(kspace.k_space_symmetry),
   k_path(kspace.k_path),
   num_nodes(kspace.num_nodes)
{
   transform_matrix = 0;
   transform_matrix += kspace.transform_matrix;
   kmesh = new libMesh::ReplicatedMesh(*(kspace.kmesh));
   kspace_comm.duplicate(kspace.kspace_comm);
}


//------------------------------------------------------------------------------//
Kspace::~Kspace()
{
  delete kmesh;
}


libMesh::MeshBase*
Kspace::get_k_mesh()
{
  return(kmesh);
}


const libMesh::MeshBase*
Kspace::get_k_mesh() const
{
  return(kmesh);
}



void
Kspace::equivalent_points(const libMesh::Point& p,
                          vector<libMesh::Point>& eq_points,
                          bool fold)
{
  Point c(p);

  // get the relative coordinates
  inverse_transform(c);


  // reduce to 1. BZ
  if (fold)
  {
    double intpart;
    c(0) = modf(c(0), &intpart);
    if (c(0) > 0.5)
      c(0) -= 1.0;
    else if (c(0) < -0.5)
      c(0) += 1.0;

    c(1) = modf(c(1), &intpart);
    if (c(1) > 0.5)
      c(1) -= 1.0;
    else if (c(1) < -0.5)
      c(1) += 1.0;

    c(2) = modf(c(2), &intpart);
    if (c(2) > 0.5)
      c(2) -= 1.0;
    else if (c(2) < -0.5)
      c(2) += 1.0;
  }

  // now c is in (-0.5, 0.5) for each direction

  eq_points.resize(0);

  switch (k_space_symmetry)
  {
    case LINEAR:
      eq_points.push_back(c);
      c(b1) *= -1;
      eq_points.push_back(c);
      break;

    case QUADRATIC:
    {
      // possible operations:
      // interchange of b1 and b2 axes + mirroring on orthogonal planes

      Point np(c);
      vector<libMesh::Point> mirrored;
      mirror(np, mirrored, {b1, b2});

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      np(b1) = c(b2);
      np(b2) = c(b1);

      mirror(np, mirrored, {b1, b2});

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      break;
    }

    case RECTANGULAR:
    {
      // possible operations:
      // mirroring on orthogonal planes

      Point np(c);
      vector<libMesh::Point> mirrored;
      mirror(np, mirrored, {b1, b2});

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      break;
    }

    case CUBIC:
    {
      // possible operations:
      // interchange of axes + mirroring on orthogonal planes

      Point np(c);
      vector<libMesh::Point> mirrored;
      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      np(0) = c(1);
      np(1) = c(0);

      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      np(0) = c(0);
      np(1) = c(2);
      np(2) = c(1);

      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      np(0) = c(2);
      np(1) = c(1);
      np(2) = c(0);

      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      break;
    }

    case TETRAGONAL:
    {
      // possible operations:
      // interchange of b1,b2 axes + mirroring on orthogonal planes

      Point np(c);
      vector<libMesh::Point> mirrored;
      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      np(b1) = c(b2);
      np(b2) = c(b1);

      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);


      break;
    }

    case ORTHORHOMBIC:
    {
      // possible operations:
      // mirroring on orthogonal planes

      eq_points.resize(0);

      Point np(c);
      vector<libMesh::Point> mirrored;
      mirror(np, mirrored);

      for (auto&& k : mirrored)
        eq_points.push_back(k);

      break;
    }

    case HEXAGONAL:
    {
      // rotations are implemented by cyclic change of
      // b1 and b2:
      // a*b1             b*b2
      // a*b2             b*(-b1 + b2)
      // a*(-b1 + b2)     b*-b1
      // a*-b1            b*-b2
      // a*-b2            b*(b1 - b2)
      // a*(b1 - b2)      b*b1
      //
      // mirroring is implemented by interchanging b1 and b2,
      // and b3 -> -b3

      // first adjust first point, because it might be outside of
      // the standard 1st BZ
      Point np(c);
      eq_points.push_back(np);

      if (fold)
      {

        //      ______
        //     /   \ /
        //    /|   |/
        //   /_\___/
        //
        //  the lower left and upper right angle are out of the
        //  standard 1st BZ, and we have to bring them back.
        //

        c(b3) = 0;
        c(b3) = np(b3);

        transform_point(np);
        double np_norm = np.norm();
        Point vec1(0);
        vec1(b1) = 1;
        transform_point(vec1);

        double prod = vec1 * np;
        double x = prod / vec1.norm();
        double y = sqrt(np_norm * np_norm - x*x);

        if ((fabs(x) > 0.5 * vec1.norm()) && (y < fabs(x) / sqrt(3)))
        {
          c(b1) -= (x > 0) ? 1 : -1;
        }
        else if ((y > fabs(x) / sqrt(3)) &&
                 (y > (vec1.norm() - fabs(x)) / sqrt(3)))
        {
          c(b2) -= (x > 0) ? 1 : -1;
        }

        eq_points[0] = c;
        np = c;
      }


      np(b1) = c(b2);
      np(b2) = c(b1);
      eq_points.push_back(np);

      np(b1) = -c(b2);
      np(b2) = c(b1) + c(b2);
      eq_points.push_back(np);

      np(b1) = -c(b1);
      np(b2) = c(b1) + c(b2);
      eq_points.push_back(np);

      np(b1) = -c(b1) - c(b2);
      np(b2) = c(b1);
      eq_points.push_back(np);

      np(b1) = -c(b1) - c(b2);
      np(b2) = c(b2);
      eq_points.push_back(np);

      np(b1) = -c(b1);
      np(b2) = -c(b2);
      eq_points.push_back(np);

      np(b1) = -c(b2);
      np(b2) = -c(b1);
      eq_points.push_back(np);

      np(b1) = c(b2);
      np(b2) = -c(b1) - c(b2);
      eq_points.push_back(np);

      np(b1) = c(b1);
      np(b2) = -c(b1) - c(b2);
      eq_points.push_back(np);

      np(b1) = c(b1) + c(b2);
      np(b2) = -c(b1);
      eq_points.push_back(np);

      np(b1) = c(b1) + c(b2);
      np(b2) = -c(b2);
      eq_points.push_back(np);

      eq_points.insert(eq_points.end(), eq_points.begin(), eq_points.end());
      for (unsigned int i = 12; i < 24; ++i)
        eq_points[i](b3) *= -1;

      break;
    }

    case FCC:
    {
      // symmetry operations:
      // rotation by 120° around \Gamma->L
      // rotation by 90° around z
      // mirror at plane \Gamma-K-L
      // mirror at plane xy
      // => 3 * 4 * 2 * 2 = 48

      Point np(c);

      // transform back to orthogonal coordinates
      transform_point(np);

      // transform to reference x,y,z
      RealTensor map(0);
      map(b1, b1) = 1; map(b1, b3) = 1;
      map(b2, b1) = 1; map(b2, b2) = 1;
      map(b3, b2) = 1; map(b3, b3) = 1;

      RealVectorValue a = map * k_basis_vector1;
      RealVectorValue b = map * k_basis_vector2;
      RealVectorValue c = map * k_basis_vector3;
      a /= a.norm();
      b /= b.norm();
      c /= c.norm();

      RealTensor Ri(a, b, c);
      RealTensor R(Ri.transpose());

      SpaceTransformation::create_star("Oh", R*np, eq_points);

      for (auto&& a : eq_points)
      {
        a = Ri*a;
        inverse_transform(a);
      }

      break;
    }

    case BCC:
      throw InitFailedException("Necessary BZ not implemented yet for unfolding");
      break;

    default:
      eq_points.push_back(Point(0));
      break;
  }

  sort(eq_points.begin(), eq_points.end());

  vector<Point>::iterator ip = unique(eq_points.begin(),
                                      eq_points.end(),
                                      compare_points);
  eq_points.resize(std::distance(eq_points.begin(), ip));

  for (auto&& a : eq_points)
    transform_point(a);

}

void
Kspace::mirror(const libMesh::Point& p,
    vector<libMesh::Point>& mirrored_points,
    const set<unsigned int>& planes)
{
  mirrored_points.resize(0);
  mirrored_points.push_back(p);

  Point np(p);
  if (planes.count(2))
  {
    np(2) *= -1;
    mirrored_points.push_back(np);
    np(2) *= -1;
  }

  if (planes.count(1))
  {
    np(1) *= -1;
    mirrored_points.push_back(np);

    if (planes.count(2))
    {
      np(2) *= -1;
      mirrored_points.push_back(np);
      np(2) *= -1;
    }
    np(1) *= -1;
  }

  if (planes.count(0))
  {
    np(0) *= -1;
    mirrored_points.push_back(np);

    if (planes.count(1))
    {
      np(1) *= -1;
      mirrored_points.push_back(np);

      if (planes.count(2))
      {
        np(2) *= -1;
        mirrored_points.push_back(np);

        np(1) *= -1;
        mirrored_points.push_back(np);
      }
    }
  }

}


void
Kspace::inverse_transform(Point& p) const
{
  Tensor1 vec;
  vec(1) = p(0);
  vec(2) = p(1);
  vec(3) = p(2);

  vec = inv(transform_matrix) * vec;

  p = Point(vec(1), vec(2), vec(3));
}

//---------------------------------------------------------------------//

void Kspace::build_k_grid()
{

  //build mesh
  kmesh = new libMesh::ReplicatedMesh(kspace_comm, k_space_dim);

  //bool full = mod_opt.get_option("full_zone", false);

  /*
   * NOTE: elements are ordered specifically, and the obtained
   * symmetry points due to axis orientation might change order.
   * In this case, which manifests as negative volumes, the node
   * order must be changed. This is implemented in some of the
   * BZ elements.
   */

  if (k_space_dim > 0)
  {

    switch (k_space_symmetry)
    {
      case LINEAR:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("X"), 1, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Edge2);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);

        break;
      }

      case QUADRATIC:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("X"), 1, 0);
        kmesh->add_point(get_symmetry_point("M"), 2, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Tri3);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);
        elem->set_node(2) = kmesh->node_ptr(2);

        break;
      }

      case RECTANGULAR:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("X"), 1, 0);
        kmesh->add_point(get_symmetry_point("S"), 2, 0);
        kmesh->add_point(get_symmetry_point("Y"), 3, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Quad4);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);
        elem->set_node(2) = kmesh->node_ptr(2);
        elem->set_node(3) = kmesh->node_ptr(3);

        break;
      }

      case CUBIC:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("X"), 1, 0);
        kmesh->add_point(get_symmetry_point("M"), 2, 0);
        kmesh->add_point(get_symmetry_point("R"), 3, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Tet4);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);
        elem->set_node(2) = kmesh->node_ptr(2);
        elem->set_node(3) = kmesh->node_ptr(3);

        break;
      }

      case TETRAGONAL:
      {
        kmesh->add_point(get_symmetry_point("M"), 0, 0);
        kmesh->add_point(get_symmetry_point("X"), 1, 0);
        kmesh->add_point(Point(0,0,0), 2, 0);
        kmesh->add_point(get_symmetry_point("A"), 3, 0);
        kmesh->add_point(get_symmetry_point("R"), 4, 0);
        kmesh->add_point(get_symmetry_point("Z"), 5, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Prism6);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);
        elem->set_node(2) = kmesh->node_ptr(2);
        elem->set_node(3) = kmesh->node_ptr(3);
        elem->set_node(4) = kmesh->node_ptr(4);
        elem->set_node(5) = kmesh->node_ptr(5);

        if (elem->volume() < 0)
        {
          elem->set_node(1) = kmesh->node_ptr(2);
          elem->set_node(2) = kmesh->node_ptr(1);
          elem->set_node(4) = kmesh->node_ptr(5);
          elem->set_node(5) = kmesh->node_ptr(4);
        }

        break;
      }

      case ORTHORHOMBIC:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("X"), 1, 0);
        kmesh->add_point(get_symmetry_point("S"), 2, 0);
        kmesh->add_point(get_symmetry_point("Y"), 3, 0);
        kmesh->add_point(get_symmetry_point("Z"), 4, 0);
        kmesh->add_point(get_symmetry_point("U"), 5, 0);
        kmesh->add_point(get_symmetry_point("R"), 6, 0);
        kmesh->add_point(get_symmetry_point("T"), 7, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Hex8);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);
        elem->set_node(2) = kmesh->node_ptr(2);
        elem->set_node(3) = kmesh->node_ptr(3);
        elem->set_node(4) = kmesh->node_ptr(4);
        elem->set_node(5) = kmesh->node_ptr(5);
        elem->set_node(6) = kmesh->node_ptr(6);
        elem->set_node(7) = kmesh->node_ptr(7);

        if (elem->volume() < 0)
        {
          elem->set_node(1) = kmesh->node_ptr(3);
          elem->set_node(3) = kmesh->node_ptr(1);
          elem->set_node(5) = kmesh->node_ptr(7);
          elem->set_node(7) = kmesh->node_ptr(5);
        }

        break;
      }

      case HEXAGONAL:
      {
        if (k_space_dim == 2)
        {
          kmesh->add_point(Point(0,0,0), 0, 0);
          kmesh->add_point(get_symmetry_point("M"), 1, 0);
          kmesh->add_point(get_symmetry_point("K"), 2, 0);

          Elem* elem = kmesh->add_elem(new libMesh::Tri3);
          elem->set_node(0) = kmesh->node_ptr(0);
          elem->set_node(1) = kmesh->node_ptr(1);
          elem->set_node(2) = kmesh->node_ptr(2);
        }
        else
        {
          kmesh->add_point(get_symmetry_point("M"), 0, 0);
          kmesh->add_point(get_symmetry_point("K"), 1, 0);
          kmesh->add_point(Point(0,0,0), 2, 0);
          kmesh->add_point(get_symmetry_point("L"), 3, 0);
          kmesh->add_point(get_symmetry_point("H"), 4, 0);
          kmesh->add_point(get_symmetry_point("A"), 5, 0);

          Elem* elem = kmesh->add_elem(new libMesh::Prism6);
          elem->set_node(0) = kmesh->node_ptr(0);
          elem->set_node(1) = kmesh->node_ptr(1);
          elem->set_node(2) = kmesh->node_ptr(2);
          elem->set_node(3) = kmesh->node_ptr(3);
          elem->set_node(4) = kmesh->node_ptr(4);
          elem->set_node(5) = kmesh->node_ptr(5);

          if (elem->volume() < 0)
          {
            elem->set_node(1) = kmesh->node_ptr(2);
            elem->set_node(2) = kmesh->node_ptr(1);
            elem->set_node(4) = kmesh->node_ptr(5);
            elem->set_node(5) = kmesh->node_ptr(4);
          }
        }

        break;
      }

      case FCC:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("K"), 1, 0);
        kmesh->add_point(get_symmetry_point("W"), 2, 0);
        kmesh->add_point(get_symmetry_point("X"), 3, 0);
        kmesh->add_point(get_symmetry_point("L"), 4, 0);
        kmesh->add_point(get_symmetry_point("U"), 5, 0);

        libMesh::TetGenMeshInterface tetgenif(*kmesh);
        tetgenif.triangulate_pointset();

        break;
      }

      case BCC:
      {
        kmesh->add_point(Point(0,0,0), 0, 0);
        kmesh->add_point(get_symmetry_point("P"), 1, 0);
        kmesh->add_point(get_symmetry_point("N"), 2, 0);
        kmesh->add_point(get_symmetry_point("H"), 3, 0);

        Elem* elem = kmesh->add_elem(new libMesh::Tet4);
        elem->set_node(0) = kmesh->node_ptr(0);
        elem->set_node(1) = kmesh->node_ptr(1);
        elem->set_node(2) = kmesh->node_ptr(2);
        elem->set_node(3) = kmesh->node_ptr(3);

        libMesh::TetGenMeshInterface tetgenif(*kmesh);
        tetgenif.triangulate_pointset();

        break;
      }

      default:
        kmesh->add_point(Point(0,0,0), 0, 0);
        break;
    }

    kmesh->prepare_for_use();
    // this is done to prevent from distributing the
    // mesh, as we want the k-mesh always serial on
    // all processes.
    kmesh->skip_partitioning(true);

    if (_mesh_order == libMesh::SECOND)
      kmesh->all_second_order();

    libMesh::MeshRefinement mr(*kmesh);

    unsigned int n = 0;
    if ((num_nodes[0] > 2) && !is_k_path())
      n = floor(log(num_nodes[0] - 1) / log(2));
    mr.uniformly_refine(n);

    //kmesh->print_info();
  }

}



//---------------------------------------------------------------------------//
void  Kspace::define_k_space(Tensor1 k_vector)
{

  Tensor1& basis1 = k_vector;

  if (basis1(1) == 1)
    transform_matrix = Tensor2(1);
  else
  {
    Tensor1 basis2;

    basis2(1) = 0;
    basis2(2) = -basis1(3);
    basis2(3) =  basis1(2);

    Tensor1 basis3 = vectorProduct(basis1, basis2);


    for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = basis2(i);
      transform_matrix(i,2) = basis3(i);
      transform_matrix(i,3) = basis1(i);

      k_basis_vector1(i-1) = basis2(i);
      k_basis_vector2(i-1) = basis3(i);
      k_basis_vector3(i-1) = basis1(i);
    }

  }

}

//-----------------------------------------------------------------------------//
void Kspace::define_k_space(Tensor1 k_vector1, Tensor1 k_vector2)
{
  // 1 -> y, 2 -> z

  Tensor1 basis1 = k_vector1;  // y
  Tensor1 basis2 = k_vector2;  // z
  Tensor1 basis3 = vectorProduct(basis1, basis2);  // x



  for (short i = 1; i < 4; i++)
  {

    transform_matrix(i,1) = basis3(i);
    transform_matrix(i,2) = basis1(i);
    transform_matrix(i,3) = basis2(i);


    k_basis_vector1(i-1) = basis3(i);
    k_basis_vector2(i-1) = basis1(i);
    k_basis_vector3(i-1) = basis2(i);

  }

  // cerr << setw(12) << transform_matrix << endl;
 
}

//----------------------------------------------------------------------------//

void Kspace::define_k_space(Tensor1 k_vector1,
                            Tensor1 k_vector2,
                            Tensor1 k_vector3)
{
 for (short i = 1; i < 4; i++)
 {
   transform_matrix(i,1) = k_vector1(i);
   transform_matrix(i,2) = k_vector2(i);
   transform_matrix(i,3) = k_vector3(i);

   k_basis_vector1(i-1) = k_vector1(i);
   k_basis_vector2(i-1) = k_vector2(i);
   k_basis_vector3(i-1) = k_vector3(i);

 }

 // cerr << setw(12) << transform_matrix << endl;

}


//----------------------------------------------------------------------------------------//
void Kspace::do_init()
{

  //if (! mod_opt.find_option("mesh_units") )
  //  throw  InitFailedException("Kspace: mesh_units must be defined");

  //double mesh_units = mod_opt.get_option("mesh_units",1.0);

  if (! mod_opt.find_option("k_space_dimension") ) 
    throw  InitFailedException("Kspace: k_space_dimension must be defined");
  
  k_space_dim = mod_opt.get_option("k_space_dimension", 0);
 
  k_path = false;
  
  mod_opt.get_option("k_max", k_max);

  if (mod_opt.find_option("k-path"))
  {
    k_path = true;
    mod_opt.get_option("number_of_nodes", num_nodes);
    if ( num_nodes.size() == 0 ) num_nodes.resize(1, 20);
  }
  else if (k_space_dim > 0)
  {
    mod_opt.get_option("number_of_nodes", num_nodes);
    if ( num_nodes.size() == 0 ) num_nodes.resize(k_space_dim, 1);
  }

  string mesh_order = mod_opt.get_option("mesh_order","first");

  if ( mesh_order == "first") _mesh_order = libMesh::FIRST;
  else if (mesh_order == "second") _mesh_order =  libMesh::SECOND;
  else throw  InitFailedException("Kspace: incorrect mesh order " + mesh_order );


  bool k_basis =  mod_opt.find_option("k1");
  //if (!k_basis && !mod_opt.find_option("r1"))
  //   throw  InitFailedException("Kspace: either k1 or r1 must be defined");
  
  std::vector<double> k_vector(3,0.0);
  

  switch (k_space_dim)
  {
    case 1:
    {
      // 1D BZ is along z

      Tensor1 vec;

      if (k_basis)
      {
        mod_opt.get_option("k1", k_vector);

        if (k_vector.size() != 3)
          throw  InitFailedException("Kspace: k1 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)
          vec(i + 1) = k_vector[i];
      }
      else
      {

        //k = 2*pi/r
        mod_opt.get_option("r1", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)  vec(i + 1) = k_vector[i];

        //vec /= (Constants::bohr_radius / mesh_units); //transform real space vector in atomic units;

        vec = (2 * M_PI * vec)/(norm(vec)*norm(vec));

      }

      define_k_space(vec);

    }
    break;

    case 2:
    {
      // 2D BZ is on y-z plane

      Tensor1 vec1;
      Tensor1 vec2;


      if (k_basis)
      {

        if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vector must be defined");

        mod_opt.get_option("k1", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: k1 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)  vec1(i + 1) = k_vector[i];

        if (! mod_opt.find_option("k2") ) throw  InitFailedException("Kspace: k2 vector must be defined");

        mod_opt.get_option("k2", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: k2 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)  vec2(i + 1) = k_vector[i];
      }
      else
      {

        Tensor1 vec1_real;
        Tensor1 vec2_real;

        mod_opt.get_option("r1", k_vector);

        if (k_vector.size() != 3)
          throw  InitFailedException("Kspace: r1 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)
          vec1_real(i + 1) = k_vector[i]; // /(Constants::bohr_radius / mesh_units);

        if (! mod_opt.find_option("r2") )
          throw  InitFailedException("Kspace: r2 vector must be defined");

        mod_opt.get_option("r2", k_vector);

        if (k_vector.size() != 3)
          throw  InitFailedException("Kspace: r2 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)
          vec2_real(i + 1) = k_vector[i]; // /(Constants::bohr_radius / mesh_units);


        Tensor1 vec3_real = vectorProduct(vec1_real, vec2_real);
        vec3_real = vec3_real/norm(vec3_real);


        double volume = (vec1_real * vectorProduct(vec2_real, vec3_real));


        if (volume == 0.0) throw  InitFailedException("Kspace: r1 and r2  vectors may be collinear ");

        vec1 = 2 * M_PI * vectorProduct(vec2_real, vec3_real) / volume;

        vec2 = 2 * M_PI * vectorProduct(vec3_real, vec1_real) / volume;

      }

      define_k_space(vec1, vec2);


    }
    break;

    case 3:
    {
      std::vector<double> k_vector;

      Tensor1 vec1;
      Tensor1 vec2;
      Tensor1 vec3;


      if (k_basis)
      {

        if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vector must be defined");

        mod_opt.get_option("k1", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: k1 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)  vec1(i + 1) = k_vector[i];

        if (! mod_opt.find_option("k2") ) throw  InitFailedException("Kspace: k2 vector must be defined");

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: k2 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)  vec2(i + 1) = k_vector[i];

        if (! mod_opt.find_option("k3") ) throw  InitFailedException("Kspace: k3 vector must be defined");

        mod_opt.get_option("k3", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: k3 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)  vec3(i + 1) = k_vector[i];

      }
      else
      {
        Tensor1 vec1_real;
        Tensor1 vec2_real;
        Tensor1 vec3_real;


        if (! mod_opt.find_option("r1") ) throw  InitFailedException("Kspace: r1 vector must be defined");

        mod_opt.get_option("r1", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)
          vec1_real(i + 1) = k_vector[i]; // /(Constants::bohr_radius / mesh_units);

        if (! mod_opt.find_option("r2") ) throw  InitFailedException("Kspace: r2 vector must be defined");

        mod_opt.get_option("r2", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: r2 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)
          vec2_real(i + 1) = k_vector[i]; // /(Constants::bohr_radius / mesh_units);

        if (! mod_opt.find_option("r3") ) throw  InitFailedException("Kspace: r3 vector must be defined");
        mod_opt.get_option("r3", k_vector);

        if (k_vector.size() != 3) throw  InitFailedException("Kspace: k3 vector size must be equal to 3");

        for (short i = 0; i < 3; i++)
          vec3_real(i + 1) = k_vector[i]; // /(Constants::bohr_radius / mesh_units);

        double volume = (vec1_real * vectorProduct(vec2_real, vec3_real));

        if (volume == 0.0) throw  InitFailedException("Kspace: r1, r2,  r3  vectors may be conplanar ");


        vec1 = 2.0 * M_PI *  vectorProduct(vec2_real, vec3_real)/volume;

        vec2 = 2.0 * M_PI *  vectorProduct(vec3_real, vec1_real)/volume;

        vec3 = 2.0 * M_PI *  vectorProduct(vec1_real, vec2_real)/volume;


      }


      define_k_space(vec1, vec2, vec3);


    }
    break;

    case 0:
      break;

    default:
      throw  InitFailedException("Kspace: k_space_dimension should be or 0 or  1 or 2 or 3");
      break;
  }

  if ((k_max(0) != 1.0) ||
      (k_max(1) != 1.0) ||
      (k_max(2) != 1.0))
  {
    ostringstream os;
    os << "scaling k space (= k_max) by (" << k_max(0) << ", "
                         << k_max(1) << ", "
                         << k_max(2) << ")";
    Messages::info(os.str());
  }



  find_k_space_symmetry();


  if (k_path)
  {
    define_k_path();
  }
  else
    build_k_grid();

  // transform the mesh to real units
  rotate_mesh();

  string filename = mod_opt.get_option("write_k_mesh", "");
  if (!filename.empty())
    libMesh::GmshIO(*kmesh).write(filename + ".msh");

}



void Kspace::find_k_space_symmetry()
{
  b1 = 0;
  b2 = 1;
  b3 = 2;

  if (k_space_dim == 0)
  {
    k_space_symmetry = GAMMA;
  }
  else if (k_space_dim == 1)
  {
    k_space_symmetry = LINEAR;
    b1 = 2;
    b2 = 0;
    b3 = 1;
  }
  else if (k_space_dim > 1)
  {
    double scalar12 = k_basis_vector1 * k_basis_vector2;
    double scalar13 = k_basis_vector1 * k_basis_vector3;
    double scalar23 = k_basis_vector2 * k_basis_vector3;

    double norm1 = k_basis_vector1.norm();
    double norm2 = k_basis_vector2.norm();
    double norm3 = k_basis_vector3.norm();

    double sp12 = scalar12/(norm1 * norm2);
    double sp13 = scalar13/(norm1 * norm3);
    double sp23 = scalar23/(norm2 * norm3);

    double angle12 = acos(sp12) * 180.0 / M_PI;
    double angle13 = acos(sp13) * 180.0 / M_PI;
    double angle23 = acos(sp23) * 180.0 / M_PI;

    // to define the symmetry of the structure angle between basis vectors are evaluated
    if (Utils::almost_equal::compare(angle12, angle13, 1e-6) &&
        Utils::almost_equal::compare(angle12, angle23, 1e-6))
    {
      //
      // all angles are the same
      //

      if (Utils::almost_equal::compare(angle12, 90, 1e-6))
      {
        //
        // orthogonal basis of BZ
        //
        if (k_space_dim == 2)
        {
          // NOTE this is not generic, assumes 2D k space on yz plane
          if (Utils::almost_equal::compare(norm2, norm3, 1e-6))
          {
            k_space_symmetry = QUADRATIC;
            b1 = 1;
            b2 = 2;
            b3 = 0;
          }
          else
          {
            k_space_symmetry = RECTANGULAR;
            b1 = 1;
            b2 = 2;
            b3 = 0;
          }
        }
        else
        {
          // Convention for tetragonal:
          //   b1,b2 -> a
          //   b3    -> c
          if (Utils::almost_equal::compare(norm1, norm2, 1e-6))
          {
            if (Utils::almost_equal::compare(norm1, norm3, 1e-6))
              k_space_symmetry = CUBIC;
            else
              k_space_symmetry = TETRAGONAL;
          }
          else if (Utils::almost_equal::compare(norm1, norm3, 1e-6))
          {
            k_space_symmetry = TETRAGONAL;
            b2 = 2;
            b3 = 1;
          }
          else if (Utils::almost_equal::compare(norm2, norm3, 1e-6))
          {
            k_space_symmetry = TETRAGONAL;
            b1 = 1;
            b2 = 2;
            b3 = 0;
          }
          else
          {
            // Convention for orthorhombic:
            //   a < b < c with a||b1, b||b2, c||b3
            k_space_symmetry = ORTHORHOMBIC;
            vector<double> len = {norm1, norm2, norm3};
            if (len[1] > len[0])
            {
              b1 = 1;
              b2 = 0;
              len[0] = norm2;
              len[1] = norm1;
            }
            if (len[2] > len[1])
            {
              b3 = b2;
              b2 = 2;
              len[2] = len[1];
              len[1] = norm3;
              len[2] = len[1];
            }
            if (len[1] > len[0])
            {
              unsigned int tmpi = b1;
              b1 = b2;
              b2 = tmpi;
              double tmpd = len[0];
              len[0] = len[1];
              len[1] = tmpd;
            }

          }
        }
      }
      else if (Utils::almost_equal::compare(sp12, -1.0/3.0, 1e-6))
      {
        k_space_symmetry = FCC;
      }
      else if (Utils::almost_equal::compare(sp12, 0.5, 1e-6))
      {
        k_space_symmetry = BCC;
      }
      else
        std::cout << "The Brillouin zone is not defined for a structure with"
        " one of the given angles between the basis vectors: "
        << angle12 << ", " << angle13 << " and " << angle23 << std::endl;

    }
    else if (Utils::almost_equal::compare(fabs(sp12), 0.5, 1e-6) &&
        Utils::almost_equal::compare(angle13, 90, 1e-5))
    {
      k_space_symmetry = HEXAGONAL;
    }
    else if (Utils::almost_equal::compare(fabs(sp13), 0.5, 1e-6) &&
        Utils::almost_equal::compare(angle12, 90, 1e-5))
    {
      b2 = 2;
      b3 = 1;
      k_space_symmetry = HEXAGONAL;
    }
    else if (Utils::almost_equal::compare(fabs(sp23), 0.5, 1e-6) &&
        Utils::almost_equal::compare(angle13, 90, 1e-5))
    {
      b1 = 1;
      b2 = 2;
      b3 = 0;
      k_space_symmetry = HEXAGONAL;
    }
    else
      std::cout << "The Brillouin zone is not defined for a structure with"
      " one of the given angles between the basis vectors: "
      << angle12 << ", " << angle13 << " and " << angle23 << std::endl;
  }

  switch (k_space_symmetry)
  {
    case LINEAR:
      degeneracy_factor = 2;
      break;

    case QUADRATIC:
      degeneracy_factor = 8;
      break;

    case RECTANGULAR:
      degeneracy_factor = 4;
      break;

    case CUBIC:
      degeneracy_factor = 32;
      break;

    case TETRAGONAL:
      degeneracy_factor = 16;
      break;

    case ORTHORHOMBIC:
      degeneracy_factor = 8;
      break;

    case HEXAGONAL:
      degeneracy_factor = 24;
      break;

    case FCC:
    case BCC:
      degeneracy_factor = 48;
      break;

    default:
      degeneracy_factor = 1;
      break;
  }

  Messages::info("Symmetry of k space: " + symmetry_names.at(k_space_symmetry));


}


libMesh::Point
Kspace::get_symmetry_point(const std::string& name) const
{
  Point p(0, 0, 0);

  // TODO it would be nicer to have the known points in some
  // static map, so that for example all known points can be
  // easily obtained, whithout duplicating code.

  if (name != "G")
  {
    switch (k_space_symmetry)
    {
      case LINEAR:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "X2")
          p(b1) = -0.5;
        break;

      case QUADRATIC:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "M")
          p(b1) = p(b2) = 0.5;
        else if (name == "X2")
          p(b2) = 0.5;
        else if (name == "M2")
        { p(b1) = -0.5; p(b2) = 0.5; }
        else if (name == "X3")
          p(b1) = -0.5;
        else if (name == "M3")
        { p(b1) = -0.5; p(b2) = -0.5; }
        else if (name == "X4")
          p(b2) = -0.5;
        else if (name == "M4")
        { p(b1) = 0.5; p(b2) = -0.5; }
        break;

      case RECTANGULAR:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "Y")
          p(b2) = 0.5;
        else if (name == "S")
          p(b1) = p(b2) = 0.5;
        else if (name == "X2")
          p(b1) = -0.5;
        else if (name == "Y2")
          p(b2) = -0.5;
        else if (name == "S2")
        { p(b1) = -0.5; p(b2) = 0.5; }
        else if (name == "S3")
        { p(b1) = -0.5; p(b2) = -0.5; }
        else if (name == "S4")
        { p(b1) = 0.5; p(b2) = -0.5; }
        break;

      case CUBIC:
        if (name == "X")
          p(0) = 0.5;
        else if (name == "M")
          p(0) = p(1) = 0.5;
        else if (name == "R")
          p(0) = p(1) = p(2) = 0.5;
        break;

      case TETRAGONAL:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "M")
          p(b1) = p(b2) = 0.5;
        else if (name == "Z")
          p(b3) = 0.5;
        else if (name == "R")
          p(b1) = p(b3) = 0.5;
        else if (name == "A")
          p(b1) = p(b2) = p(b3) = 0.5;
        break;

      case ORTHORHOMBIC:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "Y")
          p(b2) = 0.5;
        else if (name == "Z")
          p(b3) = 0.5;
        else if (name == "T")
          p(b2) = p(b3) = 0.5;
        else if (name == "U")
          p(b1) = p(b3) = 0.5;
        else if (name == "S")
          p(b1) = p(b2) = 0.5;
        else if (name == "R")
          p(b1) = p(b2) = p(b3) = 0.5;
        break;

      case HEXAGONAL:
        if (name == "A")
          p(b3) = 0.5;
        else if (name == "K")
        {
          p(b1) = 1.0/3.0;
          p(b2) = 1.0/3.0;
        }
        else if (name == "H")
        {
          p(b1) = 1.0/3.0;
          p(b2) = 1.0/3.0;
          p(b3) = 0.5;
        }
        else if (name == "M")
          p(b2) = 0.5;
        else if (name == "L")
          p(b2) = p(b3) = 0.5;
        break;

      case FCC:
        if (name == "L")
        {
          p(b1) = 0.5;
          p(b2) = 0.5;
          p(b3) = 0.5;
        }
        else if (name == "K")
        {
          p(b1) = 0.375;
          p(b2) = 0.75;
          p(b3) = 0.375;
        }
        else if (name == "U")
        {
          p(b1) = 0.25;
          p(b2) = 0.625;
          p(b3) = 0.625;
        }
        else if (name == "W")
        {
          p(b1) = 0.25;
          p(b2) = 0.75;
          p(b3) = 0.5;
        }
        else if (name == "X")
          p(b2) = p(b3) = 0.5;
        break;

      case BCC:
        if (name == "P")
        {
          p(b1) = 0.25;
          p(b2) = 0.25;
          p(b3) = 0.25;
        }
        else if (name == "N")
        {
          p(b2) = 0.5;
        }
        else if (name == "H")
        {
          p(b1) = -0.5;
          p(b2) = 0.5;
          p(b3) = -0.5;
        }
        break;

      default:
        break;
    }

    if (p == Point(0))
      throw InitFailedException("Symmetry point " + name +
          " is invalid for symmetry class " + symmetry_names[k_space_symmetry]);

    if ((k_max(0) != 1.0) ||
        (k_max(1) != 1.0) ||
        (k_max(2) != 1.0))
    {
      p(b1) *= k_max(0);
      p(b2) *= k_max(1);
      p(b3) *= k_max(2);
    }
  }

  return(p);
}

//---------------------------------------------------------------------------------------------------------------//
void Kspace::define_k_path(void)
{
  kmesh = new libMesh::ReplicatedMesh(kspace_comm, k_space_dim);


  std::string kpath = mod_opt.get_option("k_path","");
  // alternative accepted input:
  kpath = mod_opt.get_option("k-path", kpath);


  Messages::info("defining a k-path: " + kpath);
  ostringstream os;
  os <<"k-dimension: " << k_space_dim << std::endl;


  std::vector<std::string> tokens;

  Utils::tokenize(kpath, tokens, "-");

  vector<libMesh::Point> points(tokens.size());
  points[0] = get_symmetry_point(tokens[0]);

  vector<double> segments(tokens.size() - 1);
  double x = 0.0;

  for (short i = 1; i < tokens.size(); i++)
  {
    points[i] = get_symmetry_point(tokens[i]);

    libMesh::Point dp = points[i] - points[i-1];
    transform_point(dp);
    segments[i-1] = dp.norm();
    x += segments[i-1];
  }

  vector<unsigned int> n_elems(segments.size());

  os << "# segments : (";

  // we allow a special option number_of_nodes_per_segment
  if (mod_opt.find_option("number_of_nodes_per_segment") &&
      (num_nodes.size() < 2))
  {
    num_nodes.resize(0);
    num_nodes.resize(segments.size(), mod_opt.get_option("number_of_nodes_per_segment", 2));
  }

  if (num_nodes.size() > 1)
  {
    if (num_nodes.size() == segments.size())
    {
      for (short i = 0; i < segments.size(); ++i)
      {
        n_elems[i] = num_nodes[i] - 1;
        os << n_elems[i] << " ";
      }
    }
    else
    {
      throw InitFailedException("The list of number of nodes per segments "
          "in Dispersion does not match the number of k-path segments.");
    }
  }
  else
  {

    int npoints = num_nodes[0];
    for (short i = 0; i < segments.size(); ++i)
    {
      n_elems[i] = round((npoints - 1) * segments[i] / x);
      os << n_elems[i] << " ";
    }
  }

  os << ")";
  Messages::info(os.str());



  Messages m;
  m.info("Coordinates of symmetry points:");
  m.indent();
  m.info(tokens[0] + " : 0.0");


  kmesh->add_point(points[0], 0, 0);

  unsigned int id = 1;
  libMesh::Point p1(points[0]);
  x = 0.0;

  for (short i = 1; i < tokens.size(); i++)
  {
    libMesh::Point p2(points[i]);
    libMesh::Point dp = (p2 - p1);

    x = x + segments[i-1];
    ostringstream os;
    os << tokens[i] << " : " << x;
    m.info(os.str());

    dp /= n_elems[i-1];

    for (int j = 0; j < n_elems[i-1]; j++)
    {
      p1 += dp;

      kmesh->add_point(p1, id, 0);
      id++;
    }

    // we could have accumulated round off errors
    p1 = p2;
  }

  //kmesh->print_info();
}





//---------------------------------------------------------------------------------------------------------------//
void Kspace::parse_options(void)
{


}


void Kspace::transform_point(libMesh::Point& p) const
{
  Tensor1 vec1;
  vec1(1) = p(0);
  vec1(2) = p(1);
  vec1(3) = p(2);

  vec1 = transform_matrix * vec1;

  p(0) = vec1(1);
  p(1) = vec1(2);
  p(2) = vec1(3);

}

//---------------------------------------------------------------------------------------------------------------//
void Kspace::rotate_mesh(void)
{

  Tensor1 vec1;

  for (unsigned int n=0; n < kmesh->n_nodes(); n++)
  {
    const libMesh::Point p = *kmesh->node_ptr(n);
    
    vec1(1) = p(0);
    vec1(2) = p(1);
    vec1(3) = p(2);
    
    vec1 = transform_matrix * vec1;
    
    *kmesh->node_ptr(n) = libMesh::Point( vec1(1), vec1(2), vec1(3) );
    
  }

}


//---------------------------------------------------------------------------------------------------------------//
void Kspace::inv_rotate_mesh(void)
{

  Tensor1 vec1;

  Tensor2 inv_matrix = inv(transform_matrix); //.transpose();

  for (unsigned int n=0; n < kmesh->n_nodes(); n++)
  {
    const libMesh::Point p = *kmesh->node_ptr(n);
    
    vec1(1) = p(0);
    vec1(2) = p(1);
    vec1(3) = p(2);
    
    vec1 = inv_matrix * vec1;
    
    *kmesh->node_ptr(n) = libMesh::Point( vec1(1), vec1(2), vec1(3) );
    
  }

}
