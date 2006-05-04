using namespace std;
#include "SBbulkHamiltonian.h"

//=======================================================================//
void SBbulkHamiltonian::set_band_edge_energy(double energy)
{
  edge = energy/Hartree; 
}

//=======================================================================//
SBbulkHamiltonian::SBbulkHamiltonian(void)
{
  edge = 0.0;
  imass = Tensor2Sym(1);
}
//=======================================================================//
SBbulkHamiltonian::SBbulkHamiltonian(double band_edge, Tensor2Sym& imass1)
{
  edge = band_edge/Hartree;
  imass = imass1;
}
//=======================================================================//

void SBbulkHamiltonian::set_diag_mass_tensor(double m_xx, double m_yy, double m_zz)
{
  imass = Tensor2Sym(0);
  imass(1,1) = 1.0/m_xx;
  imass(2,2) = 1.0/m_yy;
  imass(3,3) = 1.0/m_zz;

 

}

//======================================================================//

void SBbulkHamiltonian::calculate_Hamiltonian_gen(void)
{
 
  //Hamiltonian is H = 1/2 * (1/m_{ij})d/dx_i d/dx_j + E0
 
  single_band_ham.constant =  edge;
  

  for (short i = 0; i < 3; i++) 
    {
      single_band_ham.linear_left[i]  = (0.0,0.0);
      single_band_ham.linear_right[i] = (0.0,0.0);
    }
  
  for (short i = 0; i < 3; i++) 
    for (short j = 0; i >= j; j++)
      {  
	single_band_ham.quad[i][j] = 0.5*imass(i+1,j+1);
	single_band_ham.quad[j][i] = 0.5*imass(i+1,j+1);
      }

   rotate_quad(single_band_ham.quad);

  
   

}

//======================================================================//
void SBbulkHamiltonian::calculate_Hamiltonian_k_par(void)
{

  MatrixElement result;
  result = single_band_ham;

 

  //------we have to change constant term
  for (short i1 = 0; i1 < 3; i1++)
    {
      result.constant += single_band_ham.linear_left[i1]  * k_vector[i1];
      result.constant += single_band_ham.linear_right[i1] * k_vector[i1];
      for (short j1 = 0; j1 < 3; j1++)
	{
	  result.constant += single_band_ham.quad[i1][j1] * k_vector[i1] * k_vector[j1]; 
	}
    }

  //------we have to change linear term
  
  for (short i1 = 0; i1 < 3; i1++)
    for (short j1 = 0; j1 < 3; j1++)
      {
	result.linear_left[i1]  += single_band_ham.quad[i1][j1] * k_vector[j1];
	result.linear_right[j1] += single_band_ham.quad[i1][j1] * k_vector[i1];
      }
  

  
  vector<MatrixElement> temp;
  temp.push_back(result); 
  Hamiltonian.push_back(temp);



  
}

//======================================================================//
//-------------------------------------------------------//
void SBbulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
  single_band_ham.constant -= el_potential/Hartree;
}


//======================================================================//
