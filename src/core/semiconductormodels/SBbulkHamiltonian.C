using namespace std;
#include "SBbulkHamiltonian.h"
#include  "ZbDDsemiconductor.h"
#include  "WzDDsemiconductor.h"
#include "getpot.h"
#include "Alloy.h"

//======================================================================//
SBbulkHamiltonian::~SBbulkHamiltonian()
{
  delete(semiconductor);
}
//=======================================================================//
SBbulkHamiltonian::SBbulkHamiltonian(const SBbulkHamiltonian& model)
{
  edge = model.edge;
  imass = model.imass;
  _filename = model._filename;
  single_band_ham = model.single_band_ham;
  semiconductor = model.semiconductor;
}

//=======================================================================//
void SBbulkHamiltonian::set_band_edge_energy(double energy)
{
  semiconductor = NULL;
  edge = energy/Hartree; 
}

//=======================================================================//
SBbulkHamiltonian::SBbulkHamiltonian(void)
{
  edge = 0.0;
  semiconductor = NULL;
  
  imass = Tensor2Sym(1);
}
//=======================================================================//
SBbulkHamiltonian::SBbulkHamiltonian(double band_edge, Tensor2Sym& imass1)
{
  edge = band_edge/Hartree;
  imass = imass1;
  semiconductor = NULL;
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

  std::vector<MatrixElement>  temp;
  temp.push_back(single_band_ham);
  Hamiltonian.push_back(temp);

  
   

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
  Hamiltonian_without_strain_pot = Hamiltonian;


  
}

//======================================================================//
//-------------------------------------------------------//
void SBbulkHamiltonian::apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential)
{
  Hamiltonian[0][0].constant = Hamiltonian_without_strain_pot[0][0].constant  -  el_potential/Hartree;
}


//======================================================================//
void SBbulkHamiltonian::read_database(const Dummy& dd)
{
  GetPot data(_filename);
  const std::string structure = data("structure", "zb");

  if (structure == "zb")
    {
      
      ZbDDsemiconductor* zbsc = new ZbDDsemiconductor();
      semiconductor = zbsc;
      
    }

  if (structure == "wz")
    {
      
      WzDDsemiconductor* wzsc = new WzDDsemiconductor();
      semiconductor =  wzsc;
    }

  semiconductor->read_database(dd);
 
  calculate_edge_and_mass();
  


}


//======================================================================//

void SBbulkHamiltonian::build_alloy(const std::string& component2,
			   const std::string& bowing_params, double content)
{
  semiconductor->build_alloy( component2, bowing_params, content);

  calculate_edge_and_mass();

}

//======================================================================//
