using namespace std;
#include "SBuserHamiltonian.h"

#include "getpot.h"

//======================================================================//
SBuserHamiltonian::~SBuserHamiltonian()
{
 
}

//=======================================================================//
SBuserHamiltonian::SBuserHamiltonian(void)
{
  edge = 0.0;
 
  
  imass = Tensor2Sym(1);
}



//======================================================================//


void SBuserHamiltonian::do_init( )
{


  EFAbulkHamiltonian::do_init();

  ModelOptions & options = get_options ();

  edge = options.get_option("Edge", 0.0);

  
  {//read (1/m) from input
    std::vector < std::vector <double>  > im;

    im.resize(3);
    
    for (short i = 0; i < 3; i++)
    {
      im[i].resize(3, 0.0);
      im[i][i] = 1.0;
    }

    options.get_option("imass", im);
    
    for (short i = 0; i < 3; i++)
      for (short j = 0; j < i; j++) 
	imass(i+1, j+1) = im[i][j];
  }

  calculate_Hamiltonian_gen();
  
}

//======================================================================//
