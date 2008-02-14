#include "Kspace.h"
#include "SimulationEnvironment.h"
using namespace std;

//------------------------------------------------------------------------------//
Kspace::Kspace()
{
  kmesh = NULL;
}


//------------------------------------------------------------------------------//
Kspace::~Kspace()
{

 
  delete kmesh;
}



//---------------------------------------------------------------------//
const Mesh& Kspace::get_k_mesh() const
{
  return(*kmesh);

}
//---------------------------------------------------------------------//

void Kspace::build_k_grid()
{


  if (k_dim > 0)
  {

    //build mesh
    kmesh = new Mesh(k_dim);



    ElemType type;
    if (integration_order == SECOND)
    {

      if (k_dim == 1)
      {
	type = EDGE3;
      }
    
      if (k_dim == 2)
      {
	type = QUAD8;
   
      }
    
      if (k_dim == 3)
      {
	type = HEX27;
      }
    }
    else 
    {
      if (k_dim == 1)
      {
	type = EDGE2;
      }
    
      if (k_dim == 2)
      {
	type = QUAD4;
   
      }
    
      if (k_dim == 3)
      {
	type = HEX8;
      }


    }

    MeshTools::Generation::build_cube (*kmesh, 
				       num_nodes[0], num_nodes[1], num_nodes[2], 
				       kmin[0], kmax[0], 
				       kmin[1], kmax[1], 
				     kmin[2], kmax[2],
				     type);



 
    rotate_mesh(kmesh, transform_matrix);
  
 
    kmesh->print_info();

 
  }

}
//---------------------------------------------------------------------------//
void  Kspace::define_k_space(Tensor1 k_vector, unsigned int n)
{


  double norm_k = norm(k_vector);

  if (wedge == HALF)
  {
    kmin[0] = 0.0; 
    kmax[0] = norm_k/2.0;
  }
  else
  {
    kmin[0] = -norm_k/2.0;  
    kmax[0] =  norm_k/2.0;
  } 
  

  num_nodes[0] = n;

  kmin[1] = 0.0; kmin[2] = 0.0;
  kmax[1] = 0.0; kmax[2] = 0.0;
  
  //TODO has to be changed!!
  Tensor1 basis1 = k_vector/norm_k;

  if (basis1(1) == 1)
    transform_matrix = Tensor2Sym(1);
  else
    {
      Tensor1 basis2;
     
      basis2(1) = 0;
      basis2(2) = -basis1(3);
      basis2(3) =  basis1(2);

      basis2 = basis2/norm(basis2);

      Tensor1 basis3 = vectorProduct(basis1, basis2);

      
      for (short i = 1; i < 4; i++)
	{
	  transform_matrix(i,1) = basis1(i);
	  transform_matrix(i,2) = basis2(i);
	  transform_matrix(i,3) = basis3(i);
	}
      
    }

 

  // k_dim = 1;

 

}

//-----------------------------------------------------------------------------//
void Kspace::define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m)
{

  double norm_k1 = norm(k_vector1);
  double norm_k2 = norm(k_vector2);

  if (wedge == QUARTER)
  {
       
    kmin[0] = 0;  kmax[0] = norm_k1/2.0; 
  
    kmin[1] = 0;  kmax[1] = norm_k2/2.0; 

  }



  if (wedge == HALF)
  {
     
    kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; 
  
    kmin[1] = 0;  kmax[1] = norm_k2/2.0; 

  }  

  if (wedge == ALL)
  {
   
    kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;

    kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;
  }

  num_nodes[0] = n; num_nodes[1] = m;

  kmin[2] = 0.0; kmax[2] = 0.0;
  
  Tensor1 basis1 = k_vector1/norm_k1;
  Tensor1 basis2 = k_vector2/norm_k2;
  Tensor1 basis3 = vectorProduct(basis1, basis2);


  for (short i = 1; i < 4; i++)
    {
      transform_matrix(i,1) = basis1(i);
      transform_matrix(i,2) = basis2(i);
      transform_matrix(i,3) = basis3(i);
    }

  //  k_dim = 2;

}

//----------------------------------------------------------------------------//

void Kspace::define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, 
				    unsigned int m, Tensor1 k_vector3, unsigned int k)
{


  double norm_k1 = norm(k_vector1);
  double norm_k2 = norm(k_vector2);
  double norm_k3 = norm(k_vector3);

 if (wedge == ALL)
 {

  
   kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;
   
   kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;

   kmin[2] = -norm_k3/2.0;  kmax[2] = norm_k3/2.0; num_nodes[2] = k;

 }

 if (wedge == HALF)
 {

   kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;

   kmin[1] = -norm_k2/2.0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;

   kmin[2] = 0;  kmax[2] = norm_k3/2.0; num_nodes[2] = k;

 }
 
 if (wedge == QUARTER)
 {

   kmin[0] = -norm_k1/2.0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;

   kmin[1] = 0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;

   kmin[2] = 0;  kmax[2] = norm_k3/2.0; num_nodes[2] = k;

 }
 

 if (wedge == EIGHTH)
 {

   kmin[0] = 0;  kmax[0] = norm_k1/2.0; num_nodes[0] = n;

   kmin[1] = 0;  kmax[1] = norm_k2/2.0; num_nodes[1] = m;

   kmin[2] = 0;  kmax[2] = norm_k3/2.0; num_nodes[2] = k;

 }
 
  
 for (short i = 1; i < 4; i++)
 {
   transform_matrix(i,1) = k_vector1(i)/norm_k1; 
   transform_matrix(i,2) = k_vector2(i)/norm_k2;
   transform_matrix(i,3) = k_vector3(i)/norm_k3;
 }



 // k_dim = 3;

}


//----------------------------------------------------------------------------------------//
void Kspace::do_init()
{
  SimulationEnvironment& si = get_environment();   

 

  double mesh_units = get_environment().get_device().get_mesh_units();
  


  const ModelOptions& mod_opt = get_options();

  if (! mod_opt.find_option("k_space_dimension") ) 
    throw  InitFailedException("Kspace: k_space_dimension must be defined");
  
  k_dim = mod_opt.get_option("k_space_dimension",0);
 
  bool k_basis =  mod_opt.get_option("k-space_basis",true);
  

  if (k_dim > 0)
  {
  
    mod_opt.get_option("number_of_nodes",num_nodes);

   

    if ( num_nodes.size() != k_dim ) 
    {
      ostringstream temp; temp << setw(4) << k_dim; 
      throw  InitFailedException("Kspace: number_of_nodes should contain " + temp.str() + " elements");
    }
  }

  {
    string wedge_type = mod_opt.get_option("wedge", "all");
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
      if (k_dim == 1)
	throw  InitFailedException("Kspace: wedge " + wedge_type + "cannot be used with 1D k-space");
    }
    else if (wedge_type == "eighth")
    {
      wedge = EIGHTH;
      degeneracy_factor = 8.0;
      if (k_dim != 3)
	throw  InitFailedException("Kspace: wedge " + wedge_type + "cannot be used only with 3D k-space");
    }
    else 
    {
      throw  InitFailedException("Kspace: unknown wedge: " + wedge_type + "\n");
    }
  }


  if (k_dim == 1)
  {

  
    

    std::vector<double> k_vector;
    Tensor1 vec;

    if (k_basis)
    {
      if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vectror must be defined"); 

      mod_opt.get_option("k1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k1 vectror size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec(i + 1) = k_vector[i];
    }
    else
    {
    
      //k = 2*pi/r
  
      if (! mod_opt.find_option("r1") ) throw  InitFailedException("Kspace: r1 vectror must be defined"); 

      mod_opt.get_option("r1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vectror size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec(i + 1) = k_vector[i];

      vec /= (Constants::bohr_radius / mesh_units); //transform real space vector in atomic units;

     
	
      vec = (2.0 * M_PI * vec)/norm(vec);
       
      

    }

    num_nodes.push_back(1);
    num_nodes.push_back(1);


    define_k_space(vec, num_nodes[0]);
    
  }
  else if (k_dim == 2)
  {
    std::vector<double> k_vector;

    Tensor1 vec1;
    Tensor1 vec2;


    if (k_basis)
    {

      if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vectror must be defined"); 

      mod_opt.get_option("k1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k1 vectror size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec1(i + 1) = k_vector[i];

      if (! mod_opt.find_option("k2") ) throw  InitFailedException("Kspace: k2 vectror must be defined"); 

      mod_opt.get_option("k2", k_vector);
      
      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k2 vectror size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec2(i + 1) = k_vector[i];
    }
    else
    {

      Tensor1 vec1_real;
      Tensor1 vec2_real;

      if (! mod_opt.find_option("r1") ) throw  InitFailedException("Kspace: r1 vectror must be defined"); 

      mod_opt.get_option("r1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vectror size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec1_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);

      if (! mod_opt.find_option("r2") ) throw  InitFailedException("Kspace: r2 vectror must be defined"); 

      mod_opt.get_option("r2", k_vector);
      
      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r2 vectror size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec2_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);


      Tensor1 vec3_real = vectorProduct(vec1_real, vec2_real);
      vec3_real = vec3_real/norm(vec3_real);


      double volume = (vec1_real * vectorProduct(vec2_real, vec3_real));
     

      if (volume == 0.0) throw  InitFailedException("Kspace: r1 and r2  vectrors may be collinear ");

      vec1 = 2.0 * M_PI* vectorProduct(vec2_real, vec3_real)/volume;

      vec2 = 2.0 * M_PI*  vectorProduct(vec3_real, vec1_real)/volume;

        

    }



    num_nodes.push_back(1);

    define_k_space(vec1, num_nodes[0],vec2, num_nodes[1] );

  }
  else if (k_dim == 3)
  {
    std::vector<double> k_vector;

    Tensor1 vec1;
    Tensor1 vec2;
    Tensor1 vec3;


    if (k_basis)
    {

      if (! mod_opt.find_option("k1") ) throw  InitFailedException("Kspace: k1 vectror must be defined");
 
      mod_opt.get_option("k1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k1 vectror size must be equal to 3");
    
      for (short i = 0; i < 3; i++)  vec1(i + 1) = k_vector[i];
      
      if (! mod_opt.find_option("k2") ) throw  InitFailedException("Kspace: k2 vectror must be defined");
      
      mod_opt.get_option("k2", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k2 vectror size must be equal to 3");
      
      for (short i = 0; i < 3; i++)  vec2(i + 1) = k_vector[i];

      if (! mod_opt.find_option("k3") ) throw  InitFailedException("Kspace: k3 vectror must be defined"); 
      mod_opt.get_option("k3", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k3 vectror size must be equal to 3");
      
      for (short i = 0; i < 3; i++)  vec3(i + 1) = k_vector[i];

      
     
     
      

    }
    else
    {
      Tensor1 vec1_real;
      Tensor1 vec2_real;
      Tensor1 vec3_real;


      if (! mod_opt.find_option("r1") ) throw  InitFailedException("Kspace: r1 vectror must be defined");
 
      mod_opt.get_option("r1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vectror size must be equal to 3");
    
      for (short i = 0; i < 3; i++)  vec1_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);
      
      if (! mod_opt.find_option("r2") ) throw  InitFailedException("Kspace: r2 vectror must be defined");
      
      mod_opt.get_option("r2", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r2 vectror size must be equal to 3");
      
      for (short i = 0; i < 3; i++)  vec2_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);

      if (! mod_opt.find_option("k3") ) throw  InitFailedException("Kspace: k3 vectror must be defined"); 
      mod_opt.get_option("k3", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k3 vectror size must be equal to 3");
      
      for (short i = 0; i < 3; i++)  vec3_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);

      
      double volume = (vec1_real * vectorProduct(vec2_real, vec3_real));
      
      if (volume == 0.0) throw  InitFailedException("Kspace: r1, r2,  r3  vectrors may be conplanar ");


      vec1 = 2.0 * M_PI *  vectorProduct(vec2_real, vec3_real)/volume;

      vec2 = 2.0 * M_PI *  vectorProduct(vec3_real, vec1_real)/volume;

      vec3 = 2.0 * M_PI *  vectorProduct(vec1_real, vec3_real)/volume;


    }
      
    
    define_k_space(vec1, num_nodes[0],vec2, num_nodes[1],vec3, num_nodes[2]);
   


  } 
  else if (k_dim == 0)
  {
    //do nothing
  } 
  else
  {
     throw  InitFailedException("Kspace: k_space_dimension should be or 0 or  1 or 2 or 3");
  }



  {
    string mesh_order  = mod_opt.get_option("mesh_order","FIRST"); ;
    if (mesh_order == "FIRST")
      integration_order =  FIRST;
    else if (mesh_order == "SECOND")
     integration_order =  SECOND;
    else throw  InitFailedException("Kspace: incorrect mesh order: " + mesh_order + "\n" );

      
  }
}
//---------------------------------------------------------------------------------------------------------------//
void Kspace::parse_options(void)
{


}

//---------------------------------------------------------------------------------------------------------------//
void Kspace::rotate_mesh (Mesh* mesh, Tensor2Gen& RotMatrix)
{

  Tensor1 vec1;

  for (unsigned int n=0; n < mesh->n_nodes(); n++)
    {
       const Point p = mesh->node(n);
  
       vec1(1) = p(0);
       vec1(2) = p(1);
       vec1(3) = p(2);

       vec1 = RotMatrix * vec1;

       

       
       
       mesh->node(n) = Point( vec1(1), vec1(2), vec1(3) );

     }

}
