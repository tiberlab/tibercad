// $Id$

#include "Kspace.h"
#include "SimulationEnvironment.h"
#include "Constants.h"
#include "Messages.h"
#include "gmsh_io.h"
#include "Utils.h"
#include <cmath>

#include <mesh_modification.h>
#include <face_tri3.h>


using namespace std;

std::vector<std::string>
Kspace::symmetry_names = {"Gamma", "linear", "quadratic", "rectangular",
                          "hexagonal", "fcc", "bcc", "cubic",
                          "tetragonal", "orthorhombic"};

//------------------------------------------------------------------------------//
Kspace::Kspace(const ModelOptions& options, const libMesh::Parallel::Communicator& comm)
 : mod_opt(options),
   _mesh_order(libMesh::FIRST)
{
  kmesh = NULL;

  kspace_comm = comm;

  transform_matrix = Tensor2Gen(1);

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
   k_space_symmetry(kspace.k_space_symmetry),
   k_path(kspace.k_path),
   num_nodes(kspace.num_nodes)
{
   transform_matrix = 0;
   transform_matrix += kspace.transform_matrix;
   kmesh = new libMesh::Mesh(*(kspace.kmesh));
   kspace_comm = kspace.kspace_comm;
}


//------------------------------------------------------------------------------//
Kspace::~Kspace()
{
  delete kmesh;
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

  if (k_space_dim > 0)
  {

    //build mesh
    kmesh = new libMesh::Mesh(kspace_comm, k_space_dim);

    libMesh::ElemType type(libMesh::EDGE2);
    if (_mesh_order == libMesh::SECOND)
    {

      if (k_space_dim == 1)
      {
        	type = libMesh::EDGE3;
      }
    
      if (k_space_dim == 2)
      {
        	type = libMesh::QUAD8;
      }
    
      if (k_space_dim == 3)
      {
        	type = libMesh::HEX27;
      }
    }
    else 
    {
      if (k_space_dim == 1)
      {
        	type = libMesh::EDGE2;
      }
    
      if (k_space_dim == 2)
      {
        	type = libMesh::QUAD4;
      }
    
      if (k_space_dim == 3)
      {
        	type = libMesh::HEX8;
      }


    }

    
    // to restrict the extension of the k-space
    libMesh::RealVectorValue k_max(1.0, 1.0, 1.0);
    mod_opt.get_option("k_max", k_max);

    // NOTE: build cube produces a mesh along x for 1D and on
    //       x-y plane for 2D, but 1D k-space is along z and 2D one
    //       is on y-z plane, so we rotate the mesh after its creation

    libMesh::MeshTools::Generation::build_cube (*kmesh,
				       num_nodes[0] - 1,
				       num_nodes[1] - 1,
				       num_nodes[2] - 1,
				       -k_max(0), k_max(0),
				       -k_max(1), k_max(1),
				       -k_max(2), k_max(2),
				       type);

    
    
    switch (k_space_dim)
    {
      case 2:
        for (unsigned int n=0; n < kmesh->n_nodes(); n++)
        {
          Point& p = kmesh->node(n);
          p(2) = p(1);
          p(1) = p(0);
          p(0) = 0.0;
        }
        // does not work for some reason:
        //MeshTools::Modification::rotate(*kmesh, 90, 90, 0);
        break;

      case 1:
        for (unsigned int n=0; n < kmesh->n_nodes(); n++)
        {
          Point& p = kmesh->node(n);
          p(2) = p(0);
          p(1) = 0.0;
          p(0) = 0.0;

        }
        //MeshTools::Modification::rotate(*kmesh, 0, 90, 90);
        break;

      default:
        break;
    }

  
  }
  else
  {
    kmesh = new libMesh::Mesh(kspace_comm, 0);
    kmesh->add_point(Point(0,0,0), 0, 0);
  }

}



//---------------------------------------------------------------------------//
void  Kspace::define_k_space(Tensor1 k_vector, unsigned int n)
{

  Tensor1& basis1 = k_vector;

  if (basis1(1) == 1)
    transform_matrix = Tensor2Sym(1);
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
void Kspace::define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m)
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

void Kspace::define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, 
				    unsigned int m, Tensor1 k_vector3, unsigned int k)
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
  

  if (mod_opt.find_option("k-path"))
  {
    k_path = true;
    mod_opt.get_option("number_of_nodes", num_nodes);
    if ( num_nodes.size() == 0 ) num_nodes.resize(1, 20);
  }


  if (k_space_dim > 0)
  {
    mod_opt.get_option("number_of_nodes", num_nodes);
    if ( num_nodes.size() == 0 ) num_nodes.resize(k_space_dim, 5);
  }

  // NOTE: for the mesh creation we always need all three indices
  // of num_nodes, so we define the missing ones as 1
  num_nodes.resize(k_space_dim);
  num_nodes.resize(3, 1);

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

      define_k_space(vec, num_nodes[0]);

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

      define_k_space(vec1, num_nodes[0],vec2, num_nodes[1] );


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


      define_k_space(vec1, num_nodes[0],vec2, num_nodes[1],vec3, num_nodes[2]);


    }
    break;

    case 0:
      break;

    default:
      throw  InitFailedException("Kspace: k_space_dimension should be or 0 or  1 or 2 or 3");
      break;
  }





  find_k_space_symmetry();


  if (k_path)
  {
    define_k_path();
  }
  else
    build_k_grid();

  //GmshIO(*kmesh).write("kspace_as_built.msh");
  // transform the mesh to real units
  rotate_mesh();

  //GmshIO(*kmesh).write("kspace.msh");

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

    double angle12 = acos(scalar12/(norm1 * norm2)) * 180.0 / M_PI;
    double angle13 = acos(scalar13/(norm1 * norm3)) * 180.0 / M_PI;
    double angle23 = acos(scalar23/(norm2 * norm3)) * 180.0 / M_PI;

    // to define the symmetry of the structure angle between basis vectors are evaluated
    if ((Utils::almost_equal::compare(angle12, 90, 1e-6)) &&
        (Utils::almost_equal::compare(angle13, 90, 1e-6)) &&
        (Utils::almost_equal::compare(angle23, 90, 1e-6)))
    {
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
          // convention: a < b < c with a||x, b||y, c||z
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
    else if (Utils::almost_equal::compare(angle12, 60, 1e-5))
    {
      k_space_symmetry = HEXAGONAL;
    }
    else if (Utils::almost_equal::compare(angle13, 60, 1e-5))
    {
      b2 = 2;
      b3 = 1;
      k_space_symmetry = HEXAGONAL;
    }
    else if (Utils::almost_equal::compare(angle23, 60, 1e-5))
    {
      b1 = 1;
      b2 = 2;
      b3 = 0;
      k_space_symmetry = HEXAGONAL;
    }

    else
      std::cout << "The Brillouin zone is not defined for a structure with"
      " one of the given angles between the basis vectors: "
      << angle12 << ", " << angle13 << " or " << angle23 << std::endl;
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

  if (name != "G")
  {
    switch (k_space_symmetry)
    {
      case LINEAR:
        if (name == "X")
          p(b1) = 0.5;
        break;

      case QUADRATIC:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "M")
          p(b1) = p(b2) = 0.5;
        break;

      case RECTANGULAR:
        if (name == "X")
          p(b1) = 0.5;
        else if (name == "Y")
          p(b2) = 0.5;
        else if (name == "S")
          p(b1) = p(b2) = 0.5;
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
      case BCC:
        break;

      default:
        break;
    }

    if (p == Point(0))
      throw InitFailedException("Symmetry point " + name +
          " is invalid for symmetry class " + symmetry_names[k_space_symmetry]);
  }

  return(p);
}

//---------------------------------------------------------------------------------------------------------------//
void Kspace::define_k_path(void)
{
  kmesh = new libMesh::Mesh(kspace_comm, k_space_dim);

  // to restrict the extension of the k-space, by an isotropic scale
  double k_max = 1.0;
  k_max = mod_opt.get_option("k_max", k_max);

  std::string kpath = mod_opt.get_option("k_path","");
  // alternative accepted input:
  kpath = mod_opt.get_option("k-path", kpath);
  int npoints = num_nodes[0];
  int nelem = max(npoints - 1, 1);


  Messages::info("(KSP) defining a k-path: " + kpath);
  ostringstream os;
  os <<"(KSP) k-dimension: "<< k_space_dim<<std::endl;
  if (k_max != 1.0)
    os <<"(KSP) scale with " << k_max << " (k_max)" << std::endl;
  os <<"(KSP) points per line "<<npoints<<std::endl; 
  Messages::info(os.str());


  std::vector<std::string> tokens;

  Utils::tokenize(kpath, tokens, "-");

  unsigned int id = 0;

  libMesh::Point p1(get_symmetry_point(tokens[0]));
  p1 *= k_max;
  kmesh->add_point(p1, id, 0);
  id++;

  for (short i = 1; i < tokens.size(); i++)
  {
    libMesh::Point p2(get_symmetry_point(tokens[i]));
    p2 *= k_max;

    libMesh::Point dp = (p2 - p1) / nelem;

    for (int j = 0; j < nelem; j++)
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

//---------------------------------------------------------------------------------------------------------------//
void Kspace::rotate_mesh(void)
{

  Tensor1 vec1;

  for (unsigned int n=0; n < kmesh->n_nodes(); n++)
  {
    const libMesh::Point p = kmesh->node(n);
    
    vec1(1) = p(0);
    vec1(2) = p(1);
    vec1(3) = p(2);
    
    vec1 = transform_matrix * vec1;
    
    kmesh->node(n) = libMesh::Point( vec1(1), vec1(2), vec1(3) );
    
  }

}


//---------------------------------------------------------------------------------------------------------------//
void Kspace::inv_rotate_mesh(void)
{

  Tensor1 vec1;

  Tensor2Gen inv_matrix = inv(transform_matrix); //.transpose();

  for (unsigned int n=0; n < kmesh->n_nodes(); n++)
  {
    const libMesh::Point p = kmesh->node(n);
    
    vec1(1) = p(0);
    vec1(2) = p(1);
    vec1(3) = p(2);
    
    vec1 = inv_matrix * vec1;
    
    kmesh->node(n) = libMesh::Point( vec1(1), vec1(2), vec1(3) );
    
  }

}
