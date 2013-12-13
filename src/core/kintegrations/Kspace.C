#include "Kspace.h"
#include "SimulationEnvironment.h"
#include "Constants.h"
#include "Messages.h"

using namespace std;

//------------------------------------------------------------------------------//
Kspace::Kspace(const ModelOptions& options)
 : mod_opt(options)
{
  kmesh = NULL;
  do_init();
}

// Copy constructor
Kspace::Kspace( const Kspace& kspace)
 : mod_opt(kspace.mod_opt),
   degeneracy_factor(kspace.degeneracy_factor),
   mesh_order(kspace.mesh_order),
   k_dim(kspace.k_dim),
   num_nodes(kspace.num_nodes)
{
   kmin[0] = kspace.kmin[0];  kmin[1] = kspace.kmin[1];  kmin[2] = kspace.kmin[2];
   kmax[0] = kspace.kmax[0];  kmax[1] = kspace.kmax[1];  kmax[2] = kspace.kmax[2];
   transform_matrix = 0;
   transform_matrix += kspace.transform_matrix;
   kmesh = new Mesh(*(kspace.kmesh));
}


//------------------------------------------------------------------------------//
Kspace::~Kspace()
{
  delete kmesh;
}


//---------------------------------------------------------------------//

void Kspace::build_k_grid()
{


  if (k_dim > 0)
  {

    //build mesh
    kmesh = new Mesh(k_dim);



    ElemType type;
    if (mesh_order == SECOND)
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

  
    rotate_mesh();
  
 
    //kmesh->print_info();

    //MeshBase::element_iterator it_k_space= kmesh->active_elements_begin();
    //const MeshBase::element_iterator it_k_end  = kmesh->active_elements_end();
 
    //for ( ; it_k_space != it_k_end ; ++it_k_space)
    // {
    //  Elem* kelem = *it_k_space;
    //  
    // cout<<"volume: " << kelem->volume()<<endl;
    //  cout<<"centroid: " << kelem->centroid()<<endl;
    //
    // }
  }
}
//---------------------------------------------------------------------------//
void  Kspace::define_k_space(Tensor1 k_vector, unsigned int n)
{


  double norm_k = norm(k_vector);

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
  

  num_nodes[0] = n-1;

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

  num_nodes[0] = n-1; num_nodes[1] = m-1;

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

  
   kmin[0] = -norm_k1;  kmax[0] = norm_k1; num_nodes[0] = n;
   
   kmin[1] = -norm_k2;  kmax[1] = norm_k2; num_nodes[1] = m;

   kmin[2] = -norm_k3;  kmax[2] = norm_k3; num_nodes[2] = k;

 }

 if (wedge == HALF)
 {

   kmin[0] = -norm_k1;  kmax[0] = norm_k1; num_nodes[0] = n;

   kmin[1] = -norm_k2;  kmax[1] = norm_k2; num_nodes[1] = m;

   kmin[2] = 0;  kmax[2] = norm_k3; num_nodes[2] = k;

 }
 
 if (wedge == QUARTER)
 {

   kmin[0] = -norm_k1;  kmax[0] = norm_k1; num_nodes[0] = n;

   kmin[1] = 0;  kmax[1] = norm_k2; num_nodes[1] = m;

   kmin[2] = 0;  kmax[2] = norm_k3; num_nodes[2] = k;

 }
 

 if (wedge == EIGHTH)
 {

   kmin[0] = 0;  kmax[0] = norm_k1; num_nodes[0] = n-1;

   kmin[1] = 0;  kmax[1] = norm_k2; num_nodes[1] = m-1;

   kmin[2] = 0;  kmax[2] = norm_k3; num_nodes[2] = k-1;

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
void Kspace::do_init() throw (InitFailedException)
{

  Messages::info("k-space creation\n");

  if (! mod_opt.find_option("mesh_units") ) 
    throw  InitFailedException("Kspace: mesh_units must be defined");

  double mesh_units = mod_opt.get_option("mesh_units",1.0);

  if (! mod_opt.find_option("k_space_dimension") ) 
    throw  InitFailedException("Kspace: k_space_dimension must be defined");
  
  k_space_dim = mod_opt.get_option("k_space_dimension",0);
 
  k_dim = k_space_dim;

  k_path = false;
  

  if (mod_opt.find_option("k-path") && k_dim > 1)
  {
  std::cout<<"(KSP) found k-path "  << std::endl;
    k_dim = 1; 
    k_path = true;
  }


  if (k_dim > 0)
  {
    //num_nodes(k_dim,0);

    mod_opt.get_option("number_of_nodes",num_nodes);
    //std::cout<<"(KSP) num nodes "<<num_nodes[0]  << std::endl;
    if ( num_nodes.size() != k_dim ) 
    {
      ostringstream temp; temp << setw(4) << k_dim; 
      throw  InitFailedException("Kspace: number_of_nodes should contain " + temp.str() + " elements");
    }

  }

  string mesh_order = mod_opt.get_option("mesh_order","first");

  if ( mesh_order == "first") mesh_order = FIRST;
  else if (mesh_order == "second") mesh_order =  SECOND;
  else throw  InitFailedException("Kspace: incorrect mesh order " + mesh_order );

  if (k_path) 
  {
    define_k_path();
    return;
  }

  {
    string def_weg;
    if (k_dim == 1) def_weg = "half";
    if (k_dim == 2) def_weg = "quarter";
    if (k_dim == 3) def_weg = "eighth";

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
      if (k_dim == 1)
	throw  InitFailedException("Kspace: wedge " + wedge_type + " cannot be used with 1D k-space");
    }
    else if (wedge_type == "eighth")
    {
      wedge = EIGHTH;
      degeneracy_factor = 8.0;
      if (k_dim != 3)
	throw  InitFailedException("Kspace: wedge " + wedge_type + " can be used only with 3D k-space");
    }
    else 
    {
      throw  InitFailedException("Kspace: unknown wedge: " + wedge_type + "\n");
    }
  }


  bool k_basis =  mod_opt.find_option("k1");
  if (!k_basis && !mod_opt.find_option("r1"))
     throw  InitFailedException("Kspace: either k1 or r1 must be defined");
  
  std::vector<double> k_vector(3,0.0);
  

  if (k_dim == 1)
  {

    Tensor1 vec;

    if (k_basis)
    {
      mod_opt.get_option("k1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k1 vector size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec(i + 1) = k_vector[i];
    }
    else
    {
    
      //k = 2*pi/r
      mod_opt.get_option("r1", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vector size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec(i + 1) = k_vector[i];

      vec /= (Constants::bohr_radius / mesh_units); //transform real space vector in atomic units;

      vec = ( M_PI * vec)/(norm(vec)*norm(vec));
       
    }

    num_nodes.push_back(0);
    num_nodes.push_back(0);


    define_k_space(vec, num_nodes[0]);
    
  }
  else if (k_dim == 2)
  {

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

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r1 vector size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec1_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);

      if (! mod_opt.find_option("r2") ) throw  InitFailedException("Kspace: r2 vector must be defined"); 

      mod_opt.get_option("r2", k_vector);
      
      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r2 vector size must be equal to 3");

      for (short i = 0; i < 3; i++)  vec2_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);


      Tensor1 vec3_real = vectorProduct(vec1_real, vec2_real);
      vec3_real = vec3_real/norm(vec3_real);


      double volume = (vec1_real * vectorProduct(vec2_real, vec3_real));
     

      if (volume == 0.0) throw  InitFailedException("Kspace: r1 and r2  vectors may be collinear ");

      vec1 = M_PI* vectorProduct(vec2_real, vec3_real)/volume;

      vec2 = M_PI* vectorProduct(vec3_real, vec1_real)/volume;
        

    }



    num_nodes.push_back(0);

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
    
      for (short i = 0; i < 3; i++)  vec1_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);
      
      if (! mod_opt.find_option("r2") ) throw  InitFailedException("Kspace: r2 vector must be defined");
      
      mod_opt.get_option("r2", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: r2 vector size must be equal to 3");
      
      for (short i = 0; i < 3; i++)  vec2_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);

      if (! mod_opt.find_option("r3") ) throw  InitFailedException("Kspace: r3 vector must be defined"); 
      mod_opt.get_option("r3", k_vector);

      if (k_vector.size() != 3) throw  InitFailedException("Kspace: k3 vector size must be equal to 3");
      
      for (short i = 0; i < 3; i++)  vec3_real(i + 1) = k_vector[i]/(Constants::bohr_radius / mesh_units);

      
      double volume = (vec1_real * vectorProduct(vec2_real, vec3_real));
      
      if (volume == 0.0) throw  InitFailedException("Kspace: r1, r2,  r3  vectors may be conplanar ");


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

  build_k_grid();
  

}
//---------------------------------------------------------------------------------------------------------------//
void Kspace::define_k_path(void)
{
  kmesh = new Mesh(k_space_dim);

  std::cout<<"(KSP) k_space_dim "<<k_space_dim<<std::endl;

  if (k_space_dim == 2)
  {    
    double G[3], M[3], X[3], X1[3];
    double *k1, *k2;
    

    double k_max = mod_opt.get_option("k_max",0.1);

    G[0]=0.0;  G[1]=0.0;  G[2]=0.0;   //( 0   0    0  ) 
    M[0]=0.0;  M[1]=0.5*k_max;  M[2]=0.5*k_max;   //( 0  1/2  1/2 )  
    X[0]=0.0;  X[1]=0.5*k_max;  X[2]=0.0*k_max;   //( 0  1/2   0  )      
    X1[0]=0.0; X1[1]=0.0; X1[2]=0.5*k_max;  //( 0   0   1/2 )        

    std::vector<std::string> tokens;
    std::string kpath = mod_opt.get_option("k_path","");
    // alternative accepted input:
    kpath = mod_opt.get_option("k-path",kpath);

    Messages::info("(KSP) defining a k-path: " + kpath);

    tokenize(kpath, tokens, "-");

    unsigned int id = 0;

    for(short i=1; i<tokens.size(); i++)
    {

      //std::cout<< "#"<<tokens[i-1]<<"-"<<tokens[i]<<"#" << endl;
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

      unsigned int npoints=num_nodes[0];

      for (unsigned int j=0; j < npoints; j++)
      {
        double b1 = k1[0]*(npoints-j)/npoints + k2[0]*j/npoints;
        double b2 = k1[1]*(npoints-j)/npoints + k2[1]*j/npoints;
        double b3 = k1[2]*(npoints-j)/npoints + k2[2]*j/npoints;

        //std::cout << b1 << " " << b2 << " " << b3 << endl;

        Point pt(b1,b2,b3);

        kmesh->add_point(pt,id,0);

        id++;
      }

    }

  }
  else if (k_space_dim == 3)
  {
    double G[3], X1[3], X2[3], X3[3], K[3], L[3];
    double *k1, *k2;
    
    double k_max = mod_opt.get_option("k_max",0.1);

    // Problems with k-points in k.p and full band!
    // Needs to define a common unit system !! 
    G[0]=0.0;    G[1]=0.0;    G[2]=0.0;  //( 0  0  0 ) 
    K[0]=0.0;    K[1]=0.5*k_max;    K[2]=0.5*k_max;  //( 0  1/2  1/2 )
    X1[0]=0.5*k_max;   X1[1]=0.0;   X1[2]=0.0;
    X2[0]=0.0;   X2[1]=0.5*k_max;   X2[2]=0.0;
    X3[0]=0.0;   X3[1]=0.0;   X3[2]=0.5*k_max;


    std::vector<std::string> tokens;
    std::string kpath = mod_opt.get_option("k-path","");

    Messages::info("(KSP) defining a k-path: " + kpath);

    tokenize(kpath, tokens, "-");

    unsigned int id = 0;

    for(short i=1; i<tokens.size(); i++)
    {

      std::cout<< "#"<<tokens[i-1]<<"-"<<tokens[i]<<"#" << endl;

       if(tokens[i-1]=="G")        k1 = G;   
       else if(tokens[i-1]=="K")   k1 = K;    
       else if(tokens[i-1]=="X1")  k1 = X1;    
       else if(tokens[i-1]=="X2")  k1 = X2;    
       else if(tokens[i-1]=="X3")  k1 = X3;   
       else if(tokens[i-1]=="L")   k1 = L;
       else                        k1 = G;  
                                              
       if(tokens[i]=="G")          k2 = G;    
       else if(tokens[i]=="K")     k2 = K;    
       else if(tokens[i]=="X1")    k2 = X1; 
       else if(tokens[i]=="X2")    k2 = X2;    
       else if(tokens[i]=="X3")    k2 = X3;   
       else if(tokens[i]=="L")     k2 = L;
       else                        k2 = G;    

      unsigned int npoints=num_nodes[0];

      for (unsigned int j=0; j < npoints; j++)
      {
        double b1 = k1[0]*(npoints-j)/npoints + k2[0]*j/npoints;
        double b2 = k1[1]*(npoints-j)/npoints + k2[1]*j/npoints;
        double b3 = k1[2]*(npoints-j)/npoints + k2[2]*j/npoints;

        //std::cout << b1 << " " << b2 << " " << b3 << endl;

        Point pt(b1,b2,b3);

        kmesh->add_point(pt,id,0);

        id++;
      }

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
    const Point p = kmesh->node(n);
    
    vec1(1) = p(0);
    vec1(2) = p(1);
    vec1(3) = p(2);
    
    vec1 = transform_matrix * vec1;
    
    kmesh->node(n) = Point( vec1(1), vec1(2), vec1(3) );
    
  }

}
//---------------------------------------------------------------------------------------------------------------//
void Kspace::inv_rotate_mesh(void)
{

  Tensor1 vec1;

  Tensor2Gen inv_matrix = transform_matrix.transpose(); 

  for (unsigned int n=0; n < kmesh->n_nodes(); n++)
  {
    const Point p = kmesh->node(n);
    
    vec1(1) = p(0);
    vec1(2) = p(1);
    vec1(3) = p(2);
    
    vec1 = inv_matrix * vec1;
    
    kmesh->node(n) = Point( vec1(1), vec1(2), vec1(3) );
    
  }

}
