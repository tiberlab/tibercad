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
   k_path(kspace.k_path),
   num_nodes(kspace.num_nodes)
{
   kmin[0] = kspace.kmin[0];  kmin[1] = kspace.kmin[1];  kmin[2] = kspace.kmin[2];
   kmax[0] = kspace.kmax[0];  kmax[1] = kspace.kmax[1];  kmax[2] = kspace.kmax[2];
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

    /* test code for MoS2
    //lattice constant for MoS2
    float a = 3.16;

    // points for triangle (1/12 of the hexagon) of the 1. Brillouin zone in y-z plane
    //G
    libMesh::Node * new_node = kmesh->add_point(libMesh::Point(0, 0, 0));
    //K
    new_node = kmesh->add_point(libMesh::Point(0, 0, 2 * M_PI/(sqrt(3) * a)));
    //M
    new_node = kmesh->add_point(libMesh::Point(0, 2 * M_PI/(3 * a), 2 * M_PI/(sqrt(3) * a)));

    //forming a triangle
    libMesh::Elem * new_elem = kmesh->add_elem(new libMesh::Tri3);

    for(int i=0; i < kmesh->n_nodes(); i++)
    {

        new_elem->set_node(i) = kmesh->node_ptr(i);

        //std::cout << kmesh->node(i) << std::endl;

    }

    kmesh->prepare_for_use();
    */


    
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
				       kmin[0] * k_max(0), kmax[0] * k_max(0),
				       kmin[1] * k_max(1), kmax[1] * k_max(1),
				       kmin[2] * k_max(2), kmax[2] * k_max(2),
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

  // k_vector is along z

  double norm_k = 0.5;

  if (wedge == HALF)
  {
    kmin[0] = 0.0; 
    kmax[0] = norm_k;
  }
  else
  {
    kmin[0] = -norm_k;  
    kmax[0] =  norm_k;
  } 
  

  kmin[1] = 0.0; kmin[2] = 0.0;
  kmax[1] = 0.0; kmax[2] = 0.0;
  
  Tensor1 basis1 = k_vector;

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

      k_basis_vector1(i-1) = basis1(i);
      k_basis_vector2(i-1) = basis2(i);
      k_basis_vector3(i-1) = basis3(i);
    }

  }

}

//-----------------------------------------------------------------------------//
void Kspace::define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m)
{
  // 1 -> y, 2 -> z

  double norm_k1 = 0.5;
  double norm_k2 = 0.5;
  //double norm_k1 = 1;
  //double norm_k2 = 1;

  if (wedge == QUARTER)
  {
       
    kmin[0] = 0;  kmax[0] = norm_k1; 
  
    kmin[1] = 0;  kmax[1] = norm_k2; 

  }



  if (wedge == HALF)
  {
     
    kmin[0] = -norm_k1;  kmax[0] = norm_k1; 
  
    kmin[1] = 0;  kmax[1] = norm_k2; 

  }  

  if (wedge == ALL)
  {
   
    kmin[0] = -norm_k1;  kmax[0] = norm_k1; num_nodes[0] = n;

    kmin[1] = -norm_k2;  kmax[1] = norm_k2; num_nodes[1] = m;
  }

  kmin[2] = 0.0; kmax[2] = 0.0;
  
  Tensor1 basis1 = k_vector1;  // y
  Tensor1 basis2 = k_vector2;  // z
  Tensor1 basis3 = vectorProduct(basis1, basis2);  // x



  for (short i = 1; i < 4; i++)
  {

    transform_matrix(i,1) = basis3(i);
    transform_matrix(i,2) = basis1(i);
    transform_matrix(i,3) = basis2(i);


    k_basis_vector1(i-1) = basis1(i);
    k_basis_vector2(i-1) = basis2(i);
    k_basis_vector3(i-1) = basis3(i);

  }

  // cerr << setw(12) << transform_matrix << endl;
 
}

//----------------------------------------------------------------------------//

void Kspace::define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, 
				    unsigned int m, Tensor1 k_vector3, unsigned int k)
{


  double norm_k1 = 0.5;
  double norm_k2 = 0.5;
  double norm_k3 = 0.5;

 if (wedge == ALL)
 {

  
   kmin[0] = -norm_k1;  kmax[0] = norm_k1;
   
   kmin[1] = -norm_k2;  kmax[1] = norm_k2;

   kmin[2] = -norm_k3;  kmax[2] = norm_k3;

 }

 if (wedge == HALF)
 {

   kmin[0] = -norm_k1;  kmax[0] = norm_k1;

   kmin[1] = -norm_k2;  kmax[1] = norm_k2;

   kmin[2] = 0;  kmax[2] = norm_k3;

 }
 
 if (wedge == QUARTER)
 {

   kmin[0] = -norm_k1;  kmax[0] = norm_k1;

   kmin[1] = 0;  kmax[1] = norm_k2;

   kmin[2] = 0;  kmax[2] = norm_k3;

 }
 

 if (wedge == EIGHTH)
 {

   kmin[0] = 0;  kmax[0] = norm_k1;

   kmin[1] = 0;  kmax[1] = norm_k2;

   kmin[2] = 0;  kmax[2] = norm_k3;

 }
 
  
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


  {
    string def_weg("all");
    switch (k_space_dim)
    {
      case 1:
        def_weg = "half";
        break;

      case 2:
        def_weg = "quarter";
        break;

      case 3:
        def_weg = "eighth";
        break;

      default:
        break;
    }

    string wedge_type = mod_opt.get_option("wedge", def_weg);
    if (wedge_type == "all" )
    {
      wedge = ALL;
      degeneracy_factor = 1.0;
    }
    else if (wedge_type == "half")
    {
      wedge = HALF;
      degeneracy_factor = 2.0;
    }
    else if (wedge_type == "quarter")
    {
      wedge = QUARTER;
      degeneracy_factor = 4.0;
      if (k_space_dim == 1)
	throw  InitFailedException("Kspace: wedge " + wedge_type + " cannot be used with 1D k-space");
    }
    else if (wedge_type == "eighth")
    {
      wedge = EIGHTH;
      degeneracy_factor = 8.0;
      if (k_space_dim != 3)
	throw  InitFailedException("Kspace: wedge " + wedge_type + " can be used only with 3D k-space");
    }
    else 
    {
      throw  InitFailedException("Kspace: unknown wedge: " + wedge_type + "\n");
    }
  }


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

    case 0:
      break;

    default:
      throw  InitFailedException("Kspace: k_space_dimension should be or 0 or  1 or 2 or 3");
      break;
  }





  define_type_of_k_space();


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


//---------------------------------------------------------------------------------------------------------------//
//calculates angular between k space basis vectors to differentiate between quadratic and hexagonal k space

void Kspace::define_type_of_k_space()
{
  if (k_space_dim == 0)
  {
    k_space_symmetry = GAMMA;
    degeneracy_factor = 1;
  }
  else if (k_space_dim == 1)
  {
    k_space_symmetry = LINEAR;
    degeneracy_factor = 2;
  }
  else if (k_space_dim > 1)
  {
    double scalar12 = k_basis_vector1* k_basis_vector2;
    double scalar13 = k_basis_vector1* k_basis_vector3;
    double scalar23 = k_basis_vector2* k_basis_vector3;

    double amountVector1 = sqrt(k_basis_vector1 * k_basis_vector1);
    double amountVector2 = sqrt(k_basis_vector2 * k_basis_vector2);
    double amountVector3 = sqrt(k_basis_vector3 * k_basis_vector3);

    double angular12 = acos(scalar12/(amountVector1 * amountVector2)) * 180.0 / M_PI;
    double angular13 = acos(scalar13/(amountVector1 * amountVector3)) * 180.0 / M_PI;
    double angular23 = acos(scalar23/(amountVector2 * amountVector3)) * 180.0 / M_PI;

    // to define the symmetry of the structure angular between basis vectors are evaluated
    if ((Utils::almost_equal::compare(angular12, 90, 1e-6)) &&
        (Utils::almost_equal::compare(angular13, 90, 1e-6)) &&
        (Utils::almost_equal::compare(angular23, 90, 1e-6)))
    {
      (amountVector1 == amountVector2 == amountVector3) ?
          k_space_symmetry = QUADRATIC : k_space_symmetry = RECTANGULAR;
    }

    //defining basis for symmetry points of the hexagonal plane
    else if (Utils::almost_equal::compare(angular12, 60, 1e-5))
    {
      identification = vector<unsigned int>(1, 3);
      k_space_symmetry = HEXAGONAL;
    }
    else if (Utils::almost_equal::compare(angular13, 60, 1e-5))
    {
      identification = vector<unsigned int>(1, 2);
      k_space_symmetry = HEXAGONAL;
    }
    else if (Utils::almost_equal::compare(angular23, 60, 1e-5))
    {
      identification = vector<unsigned int>(1, 1);
      k_space_symmetry = HEXAGONAL;
    }

    else
      std::cout << "The Brillouin zone is not defined for a structure with"
      " one of the given angulars between the basis vectors: "
      << angular12 << ", " << angular13 << " or " << angular23 << std::endl;
  }

}



//---------------------------------------------------------------------------------------------------------------//
void Kspace::define_k_path(void)
{
  kmesh = new libMesh::Mesh(kspace_comm, k_space_dim);

  // to restrict the extension of the k-space
  libMesh::RealVectorValue k_max(1.0, 1.0, 1.0);
  mod_opt.get_option("k_max", k_max);

  std::string kpath = mod_opt.get_option("k_path","");
  // alternative accepted input:
  kpath = mod_opt.get_option("k-path",kpath);
  int npoints = num_nodes[0];
  int nelem = max(npoints - 1, 1);


  Messages::info("(KSP) defining a k-path: " + kpath);
  ostringstream os;
  os <<"(KSP) k-dimension: "<< k_space_dim<<std::endl;
  os <<"(KSP) k_max: "<< k_max(0) << " " << k_max(1) << " " << k_max(2) <<std::endl;
  os <<"(KSP) points per line "<<npoints<<std::endl; 
  Messages::info(os.str());


  //note: for 2D only quadratic is implemented
  if ((k_space_dim == 2))
  {
    double G[3], M[3], X[3], X1[3];
    double *k1, *k2;


    G[0]=0.0;  G[1]=0.0;          G[2]=0.0;        //( 0   0    0  )
    M[0]=0.0;  M[1]=0.5*k_max(1); M[2]=0.5*k_max(2);  //( 0  1/2  1/2 )
    X[0]=0.0;  X[1]=0.5*k_max(1); X[2]=0.0;        //( 0  1/2   0  )
    X1[0]=0.0; X1[1]=0.0;         X1[2]=0.5*k_max(2); //( 0   0   1/2 )

    std::vector<std::string> tokens;

    tokenize(kpath, tokens, "-");

    unsigned int id = 0;

    for (short i = 1; i < tokens.size(); i++)
    {

      if(tokens[i-1]=="G")        k1 = G;
      else if(tokens[i-1]=="M")   k1 = M;
      else if(tokens[i-1]=="X")   k1 = X;
      else if(tokens[i-1]=="X'")  k1 = X1;
      else                        k1 = G;

      if(tokens[i]=="G")          k2 = G;
      else if(tokens[i]=="M")     k2 = M;
      else if(tokens[i]=="X")     k2 = X;
      else if(tokens[i]=="X'")    k2 = X1;
      else                        k2 = G;

      for (int j = (i > 1) ? 1 : 0; j < npoints; j++)
      {
        double b1 = k1[0]*(nelem-j)/nelem + k2[0]*j/nelem;
        double b2 = k1[1]*(nelem-j)/nelem + k2[1]*j/nelem;
        double b3 = k1[2]*(nelem-j)/nelem + k2[2]*j/nelem;


        libMesh::Point pt(b1,b2,b3);

        kmesh->add_point(pt,id,0);

        id++;
      }

    }

    // if there is a single point only
    if (id == 0)
    {
      if(tokens[0]=="G")          k1 = G;
      else if(tokens[0]=="M")   k1 = M;
      else if(tokens[0]=="X")   k1 = X;
      else if(tokens[0]=="X'")  k1 = X1;
      else                        k1 = G;

      kmesh->add_point(Point(k1[0], k1[1], k1[2]),id,0);
    }

  }

  else if ((k_space_dim == 3) && (k_space_symmetry == HEXAGONAL))
  {
    double G[3], M[3], K[3], A[3], L[3], H[3];
    double *k1, *k2;


    G[0]=0.0;            G[1]=0.0;               G[2]=0.0;
    M[0]=0.5;            M[1]=0.5;               M[2]=0.5;
    K[0]=0.333;          K[1]=0.333;             K[2]=0.333;
    A[0]=0.0;            A[1]=0.0;               A[2]=0.0;
    L[0]=0.5;            L[1]=0.5;               L[2]=0.5;
    H[0]=0.333;          H[1]=0.333;             H[2]=0.333;

    //basis for case 1 where k_2 and k_3 have an angular of 60 degrees:
    //   k_1  K_2  K_3
    //G   0    0    0
    //M   0    0   1/2
    //K   0   1/3  1/3
    //A  1/2   0    0
    //L  1/2   0   1/2
    //H  1/2  1/3  1/3

    std::cout << identification[0] << std::endl;
    switch (identification[0])
    {
      case 1:
        M[0] = 0;     M[1] = 0;
        K[0] = 0;
        A[0] = 0.5;
                      L[1] = 0;
        H[0] = 0.5;

        break;

      case 2:
        M[1] = 0;     M[2] = 0;
        K[1] = 0;
        A[1] = 0.5;
                      L[2] = 0;
        H[1] = 0.5;
        break;
      case 3:
        M[1] = 0;     M[2] = 0;
                      K[2] = 0;
                      A[2] = 0.5;
        L[1] = 0;
                      H[2] = 0.5;
        break;
    }
    std::vector<std::string> tokens;

    tokenize(kpath, tokens, "-");

    unsigned int id = 0;

    for (short i = 1; i < tokens.size(); i++)
    {

      if(tokens[i-1]=="G")        k1 = G;
      else if(tokens[i-1]=="M")   k1 = M;
      else if(tokens[i-1]=="K")   k1 = K;
      else if(tokens[i-1]=="A")   k1 = A;
      else if(tokens[i-1]=="L")   k1 = L;
      else if(tokens[i-1]=="H")   k1 = H;
      else                        k1 = G;

      if(tokens[i]=="G")          k2 = G;
      else if(tokens[i]=="M")     k2 = M;
      else if(tokens[i]=="K")     k2 = K;
      else if(tokens[i]=="A")     k2 = A;
      else if(tokens[i]=="L")     k2 = L;
      else if(tokens[i]=="H")     k2 = H;
      else                        k2 = G;

      for (int j = (i > 1) ? 1 : 0; j < npoints; j++)
      {
        double b1 = k1[0]*(nelem-j)/nelem + k2[0]*j/nelem;
        double b2 = k1[1]*(nelem-j)/nelem + k2[1]*j/nelem;
        double b3 = k1[2]*(nelem-j)/nelem + k2[2]*j/nelem;


        libMesh::Point pt(b1,b2,b3);

        kmesh->add_point(pt,id,0);

        id++;
      }
    }
  }

  else if ((k_space_dim == 3) && (k_space_symmetry == QUADRATIC))
  {
    double G[3], X1[3], X2[3], X3[3], M[3], M1[3], M2[3],  L[3];
    double *k1, *k2;

    G[0]=0.0;           G[1]=0.0;           G[2]=0.0;  //( 0  0  0 )
    M[0]=0.5*k_max(0);  M[1]=0.5*k_max(1);  M[2]=0.0;  //( 1/2  1/2  0 )
    M1[0]=0.5*k_max(0); M1[1]=0.0;          M1[2]=0.5*k_max(2);  //( 1/2  0  1/2 )
    M2[0]=0.0;          M2[1]=0.5*k_max(1); M2[2]=0.5*k_max(2);  //( 0  1/2  1/2 )
    X1[0]=0.5*k_max(0); X1[1]=0.0;          X1[2]=0.0;
    X2[0]=0.0;          X2[1]=0.5*k_max(1); X2[2]=0.0;
    X3[0]=0.0;          X3[1]=0.0;          X3[2]=0.5*k_max(2);
    L[0]=0.5*k_max(0);  L[1]=0.5*k_max(1);  L[2]=0.5*k_max(2);  //( 1/2  1/2  1/2 )

    std::vector<std::string> tokens;

    tokenize(kpath, tokens, "-");

    unsigned int id = 0;

    for (short i = 1; i < tokens.size(); i++)
    {


      if(tokens[i-1]=="G")        k1 = G;
      else if(tokens[i-1]=="M")   k1 = M;
      else if(tokens[i-1]=="M1")  k1 = M1;
      else if(tokens[i-1]=="M2")  k1 = M2;
      else if(tokens[i-1]=="X")   k1 = X1;
      else if(tokens[i-1]=="X1")  k1 = X1;
      else if(tokens[i-1]=="X2")  k1 = X2;
      else if(tokens[i-1]=="X3")  k1 = X3;
      else if(tokens[i-1]=="L")   k1 = L;
      else                        k1 = G;

      if(tokens[i]=="G")          k2 = G;
      else if(tokens[i]=="M")     k2 = M;
      else if(tokens[i]=="M1")    k2 = M1;
      else if(tokens[i]=="M2")    k2 = M2;
      else if(tokens[i]=="X")     k2 = X1;
      else if(tokens[i]=="X1")    k2 = X1;
      else if(tokens[i]=="X2")    k2 = X2;
      else if(tokens[i]=="X3")    k2 = X3;
      else if(tokens[i]=="L")     k2 = L;
      else                        k2 = G;


      for (int j = (i > 1) ? 1 : 0; j < npoints; j++)
      {
        double b1 = k1[0]*(nelem-j)/nelem + k2[0]*j/nelem;
        double b2 = k1[1]*(nelem-j)/nelem + k2[1]*j/nelem;
        double b3 = k1[2]*(nelem-j)/nelem + k2[2]*j/nelem;

        libMesh::Point pt(b1,b2,b3);

        kmesh->add_point(pt,id,0);

        id++;
      }

    }

    // if there is a single point only
    if (id == 0)
    {
      if(tokens[0]=="G")        k1 = G;   
      else if(tokens[0]=="M")   k1 = M; 
      else if(tokens[0]=="M1")  k1 = M1; 
      else if(tokens[0]=="M2")  k1 = M2; 
      else if(tokens[0]=="X")   k1 = X1;   
      else if(tokens[0]=="X1")  k1 = X1;    
      else if(tokens[0]=="X2")  k1 = X2;    
      else if(tokens[0]=="X3")  k1 = X3;   
      else if(tokens[0]=="L")   k1 = L;
      else                      k1 = G;  

      kmesh->add_point(Point(k1[0], k1[1], k1[2]),id,0);
    }



  }

  //kmesh->print_info();
}




//---------------------------------------------------------------------------------------------------------------//
void Kspace::tokenize(const std::string& str,
              std::vector<std::string>& tokens,
              const std::string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
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
