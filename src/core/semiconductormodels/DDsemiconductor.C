// $Id$

#include "DDsemiconductor.h"
#include "EFAbulkHamiltonian.h"
#include "Material.h"
#include "Constants.h"
#include "Database.h"


typedef std::complex<double> Complex;
extern "C"
{
  //ZHEEV( JOBZ, UPLO,  N, A,           LDA, W, WORK, LWORK, RWORK,INFO )
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& LDA, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
};



using namespace std;

//---------------------------------------------------------------------------------------------//
DDsemiconductor::~DDsemiconductor (void)
{

}
//---------------------------------------------------------------------------------------------//
DDsemiconductor::DDsemiconductor(const ModelOptions& options)
  : PhysicalModelInterface(options),
    semiconductor(NULL),
    bulk_ham(NULL),
    energy_cutoff(4.0),
    strained(false),
    k_max(1e-4)
{
}


void DDsemiconductor::prepare_submodels(void)
{
  assert(semiconductor == NULL);
  assert(bulk_ham == NULL);

  const ModelOptions& opt =  get_options();

  semiconductor = Semiconductor::create(get_material()->get_structure(), opt);
  if (semiconductor == NULL)
  {
    string msg("Cannot create semiconductor model for ");
    msg += get_material()->get_name() + " with ";
    msg += get_material()->get_structure();
    msg += " structure.";
    throw InitFailedException(msg);
  }
  add_submodel("semiconductor", semiconductor);


  ModelOptions  kp_options;
  kp_options["model"] = "kp";
  kp_options["kp_model"] = "6x6";

  bulk_ham = dynamic_cast<KPbulkHamiltonian*>(
      EFAbulkHamiltonian::create(get_material()->get_structure(), kp_options));

  if (bulk_ham == NULL)
    throw InitFailedException(string("Cannot create bulk hamiltonian for ")
        + get_material()->get_name());

  // dummy reads
  bulk_ham->get_options().get_option("kp_model", "");
  bulk_ham->get_options().get_option("model", "");

  add_submodel("bulk_ham", bulk_ham);
}

//--------------------------------------------------------------------------------------------//
void DDsemiconductor::do_init ()
{
  strain = Tensor2Sym(0);

  const ModelOptions& opt =  get_options ();

  energy_cutoff = opt.get_option("energy_cutoff", energy_cutoff);
  strained      = false;
  k_max         = opt.get_option("k_max", k_max);



}
//---------------------------------------------------------------------------------------------//
void DDsemiconductor::read_database(void)
{

  Database& db = get_database();
  db.set_section("kdotp");

  energy_cutoff  = db.get("energy_cutoff", energy_cutoff); //4eV default value
  strained       = false;
  k_max          = db.get("k_max", k_max);

}




//----------------------------------------------------------------------------------------------//

void DDsemiconductor::set_strain(const Tensor2Sym& strain_1)
{
 
  // attention: want norm() defined in tensor.h
  if (::norm( strain_1 ) > 5e-5 ) 
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

      eigvals_calculated[i] = eigvals[i]*Constants::Hartree;

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





