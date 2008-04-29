#include "Semiconductor.h"
#include "Database.h"
#include "Alloy.h"
#include "getpot.h"
#include "SimulationOptions.h"
typedef std::complex<double> Complex;
extern "C" 
{ 
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& N, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info); 
};

const double Semiconductor::Hartree = 27.2113961;


using namespace std; 

Semiconductor::Semiconductor()
{
  modelA = NULL;

  modelB = NULL;

} 

//--------------------------------------------------------------------------------------------//
void Semiconductor::do_init ()
{
  
  int verbose = SimulationOptions::verbose ();


  ModelOptions & options = get_options ();

  _consider_temperature = get_parameter("consider_temperature",  true );


/*

  if  (verbose > 0) 
    if (_consider_temperature)
      std::cout << "Semiconductor: band gap depends on temperature\n";
    else
      std::cout << "Semiconductor: band gap does not depend on temperature\n";

*/



  _temperature = get_parameter("temperature", SimulationOptions::T);





  // the temperature simulation
  string temp_simul = get_options().get_option("thermal_simulation", "");

  
  temp_interface.set_simulation(temp_simul);

/*

  if  (verbose > 0) 
    if (temp_simul != "")
      std::cout << "Semiconductor: temparature is taken from the simulation " << temp_simul << "\n";


 */ 
}





//---------------------------------------------------------------------------------------------//
KPparams   Semiconductor::calculate_kp_params (std::string kp_model )
{
  if (kp_model == "6x6")   return(calculate_6x6_kp_params());
  if (kp_model == "8x8")   return(calculate_8x8_kp_params());
}



//---------------------------------------------------------------------------------------------//
inline 
void Semiconductor::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{

  modelA = dynamic_cast<const Semiconductor* > (comp_A);

  modelB = dynamic_cast<const Semiconductor* > (comp_B);


  _xa = xa;

  do_calculate_VCA (comp_A, comp_B,  xa);


}
