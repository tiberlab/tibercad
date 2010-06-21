using namespace std;
#include "SBuserHamiltonian.h"

#include "getpot.h"

//======================================================================//
SBuserHamiltonian::~SBuserHamiltonian()
{
  
}

//=======================================================================//
SBuserHamiltonian::SBuserHamiltonian(const ModelOptions& options)
  : SBbulkHamiltonian(options)
{
  edge = 0.0;
 
  
  imass = Tensor2Sym(1);
}



//======================================================================//


void SBuserHamiltonian::do_init( )
{


  EFAbulkHamiltonian::do_init();

  edge = get_option("Edge", 0.0);

  const double mass = get_option("mass", 1.0);
 
  if (mass == 0.0) throw InitFailedException("User-defined Hamiltonian: zero mass");

  imass = (1.0 / mass) * Tensor2Sym(1.0);
  
  if (has_option("imass"))
  {//read (1/m) from input
    std::vector < std::vector <double>  > im;
    
    im.resize(3);
      
    for (short i = 0; i < 3; i++)
    {
      im[i].resize(3, 0.0);
      im[i][i] = 1.0;
    }

    get_options().get_option("imass", im);
    
    for (short i = 0; i < 3; i++)
      for (short j = 0; j < i; j++) 
	imass(i+1, j+1) = im[i][j];
  }

  calculate_Hamiltonian_gen();
  
}

//======================================================================//
