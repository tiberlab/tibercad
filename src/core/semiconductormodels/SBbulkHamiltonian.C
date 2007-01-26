using namespace std;
#include "SBbulkHamiltonian.h"

#include "getpot.h"
#include "Alloy.h"

//======================================================================//
SBbulkHamiltonian::~SBbulkHamiltonian()
{
 
}

//=======================================================================//
SBbulkHamiltonian::SBbulkHamiltonian(void)
{
 
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

  std::vector<MatrixElement>  temp;
  temp.push_back(single_band_ham);
  Hamiltonian.push_back(temp);

  Hamiltonian_without_strain_pot = Hamiltonian;
   

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
  

  

  
 
  vector<MatrixElement> temp(1,result);
  



  Hamiltonian.resize(1);

  Hamiltonian[0][0] = result;

  
  

  Hamiltonian_without_strain_pot = Hamiltonian;


 
}

//======================================================================//

void SBbulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
  Hamiltonian[0][0].constant = Hamiltonian_without_strain_pot[0][0].constant  -  el_potential/Hartree;
}


//======================================================================//

void SBbulkHamiltonian::copy_from (const PhysicalModelInterface *rhs)
{
  const SBbulkHamiltonian* mod = dynamic_cast< const SBbulkHamiltonian* > (rhs);
  
  edge = mod->edge;

  imass = mod->imass;

  single_band_ham = mod->single_band_ham;

}

//======================================================================//


void SBbulkHamiltonian::do_init( )
{

  EFAbulkHamiltonian::do_init();

  edge = 0.0;
 
  imass = Tensor2Sym(1);
  
}

//======================================================================//
