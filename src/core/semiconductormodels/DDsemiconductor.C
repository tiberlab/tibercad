
#include "DDsemiconductor.h"
typedef std::complex<double> Complex;
extern "C" 
{ 
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& N, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
};
const double DDsemiconductor::Hartree;
using namespace std;
//---------------------------------------------------------------------------------------------//
DDsemiconductor::DDsemiconductor()
{
  strain = Tensor2Sym(0);
  energy_cutoff=1.0; //1eV default value
  strained = false;
  k_max = 1e-4;
}


//--------------------------------------------------------------------------------------------//

DDsemiconductor::DDsemiconductor(const  Tensor2Sym& strain_1, const double energy_cutoff_1)
{

  if (norm( strain_1 ) > 1e-5 ) 
    {
      strain = strain_1;
      strained = true;
    }
  else
    {
      strain = Tensor2Sym(0);
      strained = false;
    }
  
  energy_cutoff = energy_cutoff_1;
  k_max = 1e-4;
}
//----------------------------------------------------------------------------------------------//
void DDsemiconductor::set_strain(const Tensor2Sym& strain_1)
{
  
 
  if (norm( strain_1 ) > 1e-5 ) 
    {
      strain = strain_1;
      strained = true;
    }
  else
    {
      strain = Tensor2Sym(0);
      strained = false;
    }
}
//---------------------------------------------------------------------------------------------//
const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_conduction_band_energy_mass(void) const
{
  // const std::vector<DDsemiconductor::band_extremum>&  result;

  // result = &conduction_band;

  return(conduction_band);

}
//---------------------------------------------------------------------------------------------//
const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_valence_band_energy_mass(void) const
{
  
  return(valence_band);

  
  

}
//---------------------------------------------------------------------------------------------//

vector< vector<double> > DDsemiconductor::calculate_vb_bulk_states(const vector<Tensor1>& k_vector)
{
 

  vector< vector<double> > result;

  KPbulkHamiltonian  bulk("6x6");
  bulk.strainM = strain;

  //cerr << setw(16) << bulk.strainM ;
 
  KPbulkHamiltonian::KPparams params_kp ;
 
  params_kp =  calculate_6x6_kp_params();
  
  
  
  bulk.kpVVtermSymmetric = true;


  bulk.set_parameters( params_kp );

  Tensor2Gen RotMatrix;

  RotMatrix =  Tensor2Gen(1);

  bulk.set_rotation_matrix( RotMatrix );

  bulk.calculate_Hamiltonian_gen();

  

  double kvec[3];


  std::complex<double> ham6x6matrix[6*6];

  unsigned int N = k_vector.size();

  vector <double> eigvals_calculated(6);
  
  for (short i1 = 0; i1 < N; i1++ )

    {

      //      kvec[0] =  k_vector[i1](1);  kvec[1] =  k_vector[i1](2);  kvec[2] =  k_vector[i1](3);


      bulk.set_k_vector(k_vector[i1]);

      bulk.calculate_Hamiltonian_k_par();

      

      std::vector<std::vector<KPbulkHamiltonian::MatrixElement> >&    Ham1 =  bulk.get_Hamiltonian() ;

      

      for (short i = 0; i < 6; i++)
	{
	  for (short j = 0; j < 6; j++)
	    {
	     
	      ham6x6matrix[i + j*6] = Ham1[i][j].constant;
	     
	    }
	}

    
      char jobs = 'N';
      char UPLO = 'U'; 
      int  N = 6;
      double eigvals[6];
      std::complex<double> WORK[11];
      int LWORK = 11;
      double RWORK[16];
      int info;
      
      


      zheev_(jobs, UPLO, N, ham6x6matrix, N, eigvals, WORK, LWORK, RWORK, info); 
      if (info !=0 ) exit(1);

     
      for (short i = 0; i < 6; i++)
	{
	  
	  eigvals_calculated[i] = eigvals[i]*Hartree;
	 
	}

     

      result.push_back(eigvals_calculated);
      

    }


  return(result);
}

//--------------------------------------------------------//


vector<vector<double> > DDsemiconductor::get_valence_kp_dispersion(Tensor1 k_i, Tensor1 k_f, unsigned int number_of_points)
{
  if (number_of_points < 2) number_of_points = 2;

  Tensor1 dk = (k_f - k_i)/(number_of_points - 1);

  vector<Tensor1> k_points(number_of_points);

  
  for (unsigned int point = 0 ; point < number_of_points; point++)
    {
      k_points[point] = k_i +  point * dk;
     
    }


  vector<vector<double> > result = calculate_vb_bulk_states(k_points);

 
  

  return(result);

}


//---------------------------------------------------------------------------------------------//
DDsemiconductor::~DDsemiconductor (void)
{

}
//---------------------------------------------------------------------------------------------//
