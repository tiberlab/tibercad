// $Id$

#include "tibercad/physics/semiconductormodels/EFAbulkHamiltonian.h"
#include "tibercad/physics/Material.h" 
#include "tibercad/math/Tensor1.h" 

using namespace std; 

//------------------------------------------------------------//
void EFAbulkHamiltonian::do_init()
{

  std::vector <double> k(3, 0.0);
  ModelOptions & options = get_options ();
  options.get_option ("k", k);

  k_vector[0] =  k[0];
  k_vector[1] =  k[1];
  k_vector[2] =  k[2];

  _degeneracy = 1;

  set_rotation_matrix();

  get_option("particle","");
  get_option("spurious","");
  get_option("temperature_scaling","");
 
}



//-------------------------------------------------------------//
EFAbulkHamiltonian::EFAbulkHamiltonian(const ModelOptions& options)
 : PhysicalModel(options)
{
 

}

EFAbulkHamiltonian* EFAbulkHamiltonian::create (const Material* mat,  const ModelOptions& options)
{

  if (! (options.find_option("model")) )
  {
    throw InitFailedException("A model must be defined for the EFA bulk hamiltonian.");
  }

  const std::string&  model_name = options.get_option("model", "8x8");
  const std::string&  particle = options.get_option("particle", "");

  std::string model;
  std::string structure = mat->get_structure();

  if (model_name == "single_band")
  {
    if (particle == "el")
      model ="quantum_cond_band_" + structure;
    else if (particle == "hl")
      model = "quantum_user";
  }
  else if (model_name == "conduction_band")
    model ="quantum_cond_band_" + structure;
  else if (model_name == "valence_band")
    model = "quantum_user";
  else
    model = "quantum_kp";

  return PhysicalModel::create<EFAbulkHamiltonian>(model, mat, options);

}
 
//------------------------------------------------------------//

void EFAbulkHamiltonian::set_k_vector (const double k_vector_in[3])
{
 k_vector[0] =  k_vector_in[0];
 k_vector[1] =  k_vector_in[1];
 k_vector[2] =  k_vector_in[2];
}


void EFAbulkHamiltonian::set_k_vector (const Tensor1& k_vector_in)
{
  k_vector[0] =  k_vector_in(1);
  k_vector[1] =  k_vector_in(2);
  k_vector[2] =  k_vector_in(3);


}


//---------------------------------------------------------------//


std::vector< std::vector<EFAbulkHamiltonian::MatrixElement > >& EFAbulkHamiltonian::get_Hamiltonian(void)
{
  return (Hamiltonian);
}


void EFAbulkHamiltonian::set_rotation_matrix()

{

  const Material* mat =	get_material();
  const libMesh::RealTensor& rotmatrix = mat->get_rotation_matrix();
  
  for (short i = 0; i < 3; i++)
     for (short j = 0; j < 3; j++)
       rot_matrix[i][j] = rotmatrix(i,j);

}


//--------------------------------------------------------//
void EFAbulkHamiltonian::rotate_linear(std::complex<double> *vector)
{
  complex<double> vec1[3];

  for (short i = 0 ; i < 3; i++) vec1[i] = vector[i];

  for (short i = 0 ; i < 3; i++)
  {
    vector[i] = 0.0;
    for (short j = 0 ; j < 3; j++)
      vector[i] += vec1[j] * rot_matrix[i][j];
  }
	
}

//-------------------------------------------------------//
void EFAbulkHamiltonian::rotate_quad(std::complex<double> matrix[][3])
{

  


  complex<double> mat1[3][3];
  for (short i = 0 ; i < 3; i++)
    for (short j = 0 ; j < 3; j++) 
      mat1[i][j] = matrix[i][j];


  for (short i = 0 ; i < 3; i++)
    for (short j = 0 ; j < 3; j++)
      matrix[i][j] =  complex<double>(0.0, 0.0);


  for (short i = 0 ; i < 3; i++)
    for (short j = 0 ; j < 3; j++)
      for (short i1 = 0 ; i1 < 3; i1++)
	for (short j1= 0 ; j1 < 3; j1++)
	  {
	    matrix[i][j] += rot_matrix[i][i1] * rot_matrix[j][j1] * mat1[i1][j1];
	  }
 
  
 
     
}
//--------------------------------------------------------//

double
EFAbulkHamiltonian::get_degeneracy(void)
{
  return _degeneracy;
}



//---------------------------------------------------------// 
const std::map<short, short>&  EFAbulkHamiltonian::get_kp_bands_map(void) const
{
  return(kp_bands_map);
}
