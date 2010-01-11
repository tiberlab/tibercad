// $Id$

#include "EFAbulkHamiltonian.h"
#include "RotatedCrystal.h"
#include "Material.h" 
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


  set_rotation_matrix();
 
}



//-------------------------------------------------------------//
EFAbulkHamiltonian::EFAbulkHamiltonian(const ModelOptions& options)
 : PhysicalModelInterface(options)
{
 

}
 
//------------------------------------------------------------//

void EFAbulkHamiltonian::set_k_vector (const double k_vector_in[3])
{
 k_vector[0] =  k_vector_in[0];
 k_vector[1] =  k_vector_in[1];
 k_vector[2] =  k_vector_in[2];
}


void EFAbulkHamiltonian::set_k_vector (Tensor1 k_vector_in)
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
  
  const RotatedCrystal& cr = mat->get_rotated_crystal ();

  Tensor2Gen rotmatrix = cr.RotMatrix;

  for (short i = 0; i < 3; i++)
     for (short j = 0; j < 3; j++)
       rot_matrix[i][j] = rotmatrix(i+1,j+1);

}


//--------------------------------------------------------//
void EFAbulkHamiltonian::rotate_linear(std::complex<double> *vector)
{
  complex<double> vec1[3];

  for (short i = 0 ; i < 3; i++) vec1[i] = vector[i];

  for (short i = 0 ; i < 3; i++)  vector[i] = complex<double>(0.0,0.0);

  for (short i = 0 ; i < 3; i++)
    for (short j = 0 ; j < 3; j++)
      vector[i] += vec1[j] * rot_matrix[i][j];
	
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





//---------------------------------------------------------// 
const std::map<short, short>&  EFAbulkHamiltonian::get_kp_bands_map(void) const
{
  return(kp_bands_map);
}
