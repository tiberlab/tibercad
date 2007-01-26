#include "Semiconductor.h"
#include "Database.h"
#include "Alloy.h"
#include "getpot.h"
typedef std::complex<double> Complex;
extern "C" 
{ 
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& N, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info); 
};

const double Semiconductor::Hartree = 27.2113961;


using namespace std; 

 

//--------------------------------------------------------------------------------------------//
void Semiconductor::do_init ()
{
 
}





//---------------------------------------------------------------------------------------------//
KPparams   Semiconductor::calculate_kp_params (std::string kp_model )
{
  if (kp_model == "6x6")   return(calculate_6x6_kp_params());
  if (kp_model == "8x8")   return(calculate_8x8_kp_params());
}



//---------------------------------------------------------------------------------------------//
