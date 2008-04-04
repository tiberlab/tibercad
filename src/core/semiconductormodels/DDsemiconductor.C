#include "DDsemiconductor.h"
#include "EFAbulkHamiltonian.h"
#include "Database.h"
#include "Alloy.h"
#include "getpot.h"
typedef std::complex<double> Complex;
extern "C" 
{ 
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& N, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info); 
};

const double DDsemiconductor::Hartree = 27.2113961;


using namespace std; 

//---------------------------------------------------------------------------------------------//
DDsemiconductor::~DDsemiconductor (void)
{
  PhysicalModelInterface::destroy(semiconductor);

  PhysicalModelInterface::destroy(bulk_ham);
}
//---------------------------------------------------------------------------------------------//
DDsemiconductor::DDsemiconductor (void)
{
  semiconductor = NULL;

  bulk_ham = NULL;
}
 

//--------------------------------------------------------------------------------------------//
void DDsemiconductor::do_init ()
{
  strain = Tensor2Sym(0);
 
  const ModelOptions& opt =  get_options ();
  
  energy_cutoff = opt.get_option("energy_cutoff", 4.0);
  strained      = false;
  k_max         = opt.get_option("k_max", 1e-3);

  PhysicalModelInterface::destroy(semiconductor);
  
  semiconductor = Semiconductor::create(get_material()->get_structure(), opt);
  if (semiconductor == NULL)
  {
    string msg("DDsemiconductor: cannot create model for material with ");
    msg += get_material()->get_structure();
    msg += " structure.";
    throw InitFailedException(msg);
  }

  semiconductor->set_material(get_material());

  semiconductor->init();
  
  PhysicalModelInterface::destroy(bulk_ham);

  ModelOptions  kp_options; 
  kp_options["model"] = "kp";
  kp_options["kp_model"] = "6x6";

  bulk_ham = dynamic_cast<KPbulkHamiltonian*>(
      EFAbulkHamiltonian::create(get_material()->get_structure(), kp_options));

  if (bulk_ham == NULL)
    throw InitFailedException("DDsemiconductor::do_init ()     bulk_ham == NULL ");

  bulk_ham->set_material(get_material());
   
  bulk_ham->init();
  

}
//---------------------------------------------------------------------------------------------//
void DDsemiconductor::read_database(void)
{

  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  energy_cutoff  = data("energy_cutoff", 4.0); //4eV default value
  strained       = false;
  k_max          = data("k_max", 1e-5);

}
//---------------------------------------------------------------------------------------------//
void DDsemiconductor::copy_from (const PhysicalModelInterface *rhs)
{
    const DDsemiconductor* mod = dynamic_cast<const DDsemiconductor*>(rhs);
  
    strain = mod->strain;
    energy_cutoff = mod->energy_cutoff;
    
    k_max = mod->k_max;

}

//--------------------------------------------------------------------------------------------//
void DDsemiconductor::calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa)
{
 
  const DDsemiconductor* modA = dynamic_cast<const DDsemiconductor*>(comp_A);
  const DDsemiconductor* modB = dynamic_cast<const DDsemiconductor*>(comp_B);
 
 
  
  energy_cutoff = alloy(modA->energy_cutoff,modB->energy_cutoff,xa);
   
  k_max = alloy(modA->k_max, modB->k_max, xa);

  semiconductor -> build_alloy(modA->semiconductor, modB->semiconductor, xa);
  bulk_ham -> build_alloy(modA->bulk_ham, modB->bulk_ham, xa);
  
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

 

  bulk_ham->apply_strain_and_potential(strain, 0.0);



}

 

//---------------------------------------------------------------------------------------------//
const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_conduction_band_energy_mass(void) const
{

   return(conduction_band);

} 
//---------------------------------------------------------------------------------------------//
const std::vector<DDsemiconductor::band_extremum>& DDsemiconductor::get_valence_band_energy_mass(void) const
{

  cerr << "dd " << valence_band.size() << "\n";
  
  return(valence_band);

}

//---------------------------------------------------------------------------------------------//

vector< vector<double> > DDsemiconductor::calculate_vb_bulk_states(const vector<Tensor1>& k_vector)
{
 
  vector< vector<double> > result;  
 
  double kvec[3];


  std::complex<double> ham6x6matrix[6*6];
   
  unsigned int N = k_vector.size();
   
  vector <double> eigvals_calculated(6);
   
  for (short i1 = 0; i1 < N; i1++ )

  {
     
  

    bulk_ham->set_k_vector(k_vector[i1]);
     
    bulk_ham->calculate_Hamiltonian_k_par();

    if (strained)  bulk_ham->apply_strain_and_potential(strain, 0.0); 

    std::vector<std::vector<KPbulkHamiltonian::MatrixElement> >&    Ham1 =  bulk_ham->get_Hamiltonian() ;
     
     
     
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





