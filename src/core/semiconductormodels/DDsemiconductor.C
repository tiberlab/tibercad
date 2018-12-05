// $Id$

#include "DDsemiconductor.h"
#include "EFAbulkHamiltonian.h"
#include "RotatedCrystal.h"
#include "Material.h"
#include "Constants.h"
#include "Database.h"
#include "Messages.h"
#include "SolverException.h"
#include "RuntimeException.h"

#include "tensor_value.h"
#include "dense_matrix.h"


typedef std::complex<double> Complex;

namespace {
  enum KDirection
  {
    X, Y, Z, XY, XZ, YZ
  };



  void rotate(libMesh::DenseVector<Complex>& vec, const RotatedCrystal& cr)
  {
    // NOTE This is not tested well
    libMesh::DenseVector<Complex> tmp(vec);

    const Tensor2Gen& mat = cr.RotMatrix;

    double a, b, c; // euler angles
    cr.get_euler_angles(a, b, c);
    a /= 2.0;
    b /= 2.0;
    c /= 2.0;

    double cosb = cos(b);
    double sinb = sin(b);
    Complex eacp = Complex(cos(a + c), sin(a + c));
    Complex eacm = Complex(cos(a - c), sin(a - c));
    Complex uu(cosb * eacp);
    Complex ud(sinb * conj(eacm));
    Complex du(-sinb * eacm);
    Complex dd(cosb * conj(eacp));

    for (int i = 0; i < 3; ++i)
    {
      int I = i + 1;
      //vec(i) = uu * (mat(I,1) * tmp(0) + mat(I,2) * tmp(1) + mat(I,3) * tmp(2)) +
      //    ud * (mat(I,1) * tmp(3) + mat(I,2) * tmp(4) + mat(I,3) * tmp(5));
      //vec(i+3) = dd * (mat(I,1) * tmp(3) + mat(I,2) * tmp(4) + mat(I,3) * tmp(5)) +
      //    du * (mat(I,1) * tmp(0) + mat(I,2) * tmp(1) + mat(I,3) * tmp(2));

      //vec(i) = conj(uu) * (mat(1,I) * tmp(0) + mat(2,I) * tmp(1) + mat(3,I) * tmp(2)) +
      //    conj(du) * (mat(1,I) * tmp(3) + mat(2,I) * tmp(4) + mat(3,I) * tmp(5));
      //vec(i+3) = conj(dd) * (mat(1,I) * tmp(3) + mat(2,I) * tmp(4) + mat(3,I) * tmp(5)) +
      //    conj(ud) * (mat(1,I) * tmp(0) + mat(2,I) * tmp(1) + mat(3,I) * tmp(2));

      vec(i) = (uu) * (mat(1,I) * tmp(0) + mat(2,I) * tmp(1) + mat(3,I) * tmp(2)) +
          (ud) * (mat(1,I) * tmp(3) + mat(2,I) * tmp(4) + mat(3,I) * tmp(5));
      vec(i+3) = (dd) * (mat(1,I) * tmp(3) + mat(2,I) * tmp(4) + mat(3,I) * tmp(5)) +
          (du) * (mat(1,I) * tmp(0) + mat(2,I) * tmp(1) + mat(3,I) * tmp(2));
    }
  }


  void rotate(libMesh::DenseVector<Complex>& vec, KDirection dir)
  {
    libMesh::DenseVector<Complex> tmp(vec);
    switch (dir)
    {
      case X:
        vec(0) = -(tmp(2) + tmp(5));
        vec(1) = tmp(1) + tmp(4);
        vec(2) = tmp(0) + tmp(3);
        vec(3) = tmp(2) - tmp(5);
        vec(4) = tmp(4) - tmp(1);
        vec(5) = tmp(3) - tmp(0);
        vec.scale(M_SQRT1_2);
        break;

      case Y:
        vec(0) = -Complex(1,1)*tmp(2) - Complex(1,-1)*tmp(5);
        vec(1) = -Complex(1,1)*tmp(0) - Complex(1,-1)*tmp(3);
        vec(2) = Complex(1,1)*tmp(1) + Complex(1,-1)*tmp(4);
        vec(3) = Complex(1,1)*tmp(2) - Complex(1,-1)*tmp(5);
        vec(4) = Complex(1,1)*tmp(0) - Complex(1,-1)*tmp(3);
        vec(5) = -Complex(1,1)*tmp(1) + Complex(1,-1)*tmp(4);
        vec.scale(0.5);
        break;

      default:
        break;

      case XY:
      {
        const Complex ep(0.923879532511287, 0.382683432365090);
        const Complex em(0.923879532511287, -0.382683432365090);
        vec(0) = -M_SQRT2 * (ep * tmp(2) + em * tmp(5));
        vec(1) = ep * (tmp(1) - tmp(0)) + em * (tmp(4) - tmp(3));
        vec(2) = ep * (tmp(1) + tmp(0)) + em * (tmp(4) + tmp(3));
        vec(3) = M_SQRT2 * (ep * tmp(2) - em * tmp(5));
        vec(4) = ep * (tmp(0) - tmp(1)) + em * (tmp(4) - tmp(3));
        vec(5) = -ep * (tmp(0) + tmp(1)) + em * (tmp(4) + tmp(3));
        vec.scale(0.5);
        break;
      }

      case XZ:
      {
        const double sinpi_8 = 0.382683432365090;
        const double cospi_8 = 0.923879532511287;
        vec(0) = cospi_8 * (tmp(0) - tmp(2)) + sinpi_8 * (tmp(3) - tmp(5));
        vec(1) = cospi_8 * tmp(1) + sinpi_8 * tmp(4);
        vec(2) = cospi_8 * (tmp(0) + tmp(2)) + sinpi_8 * (tmp(3) + tmp(5));
        vec(3) = -sinpi_8 * (tmp(0) - tmp(2)) + cospi_8 * (tmp(3) - tmp(5));
        vec(4) = -sinpi_8 * tmp(1) + cospi_8 * tmp(4);
        vec(5) = -sinpi_8 * (tmp(0) + tmp(2)) + cospi_8 * (tmp(3) + tmp(5));
        vec.scale(M_SQRT1_2);
        break;
      }

      case YZ:
      {
        const double sinpi_8 = 0.382683432365090;
        const double cospi_8 = 0.923879532511287;
        const Complex uu = cospi_8 * Complex(1,1);
        const Complex ud = sinpi_8 * Complex(1,-1);
        const Complex du = -sinpi_8 * Complex(1,1);
        const Complex dd = cospi_8 * Complex(1,-1);
        vec(0) = uu * (tmp(1) - tmp(2)) + ud * (tmp(4) - tmp(5));
        vec(1) = -M_SQRT2 * (uu * tmp(0) + ud * tmp(3));
        vec(2) = uu * (tmp(1) + tmp(2)) + ud * (tmp(4) + tmp(5));
        vec(3) = du * (tmp(1) - tmp(2)) + dd * (tmp(4) - tmp(5));
        vec(4) = -M_SQRT2 * (du * tmp(0) + dd * tmp(3));
        vec(5) = du * (tmp(1) + tmp(2)) + dd * (tmp(4) + tmp(5));
        vec.scale(0.5);
        break;
      }
    }
  }
}




extern "C"
{
  void zheev_(char& jobs, char& UPLO, int& N, Complex* ham6x6matrix, int& LDA, double* eigvals,  Complex* WORK, int& LWORK, double* RWORK, int& info);
}



using namespace std;

//---------------------------------------------------------------------------------------------//
DDsemiconductor::~DDsemiconductor (void)
{

}
//---------------------------------------------------------------------------------------------//
DDsemiconductor::DDsemiconductor(const ModelOptions& options)
  : PhysicalModel(options),
    semiconductor(NULL),
    bulk_ham(NULL),
    energy_cutoff(10.0),
    strained(false),
    k_max(1e-4)    //k_max(0.0529177e-4)
{
}


DDsemiconductor* DDsemiconductor::create(const Material* mat,
    const ModelOptions& options)
{
  std::string structure = mat->get_structure();
  return PhysicalModel::create<DDsemiconductor>("DDsemicond_" + structure, mat, options);
}

void DDsemiconductor::prepare_submodels(void)
{
  assert(semiconductor == NULL);
  assert(bulk_ham == NULL);

  //create semiconductor for electron parameters
  ModelOptions  kp_options(get_options());
  kp_options["model"] = "single_band";

  semiconductor = Semiconductor::create(get_material(), kp_options);
  if (semiconductor == NULL)
  {
    string msg("Cannot create semiconductor model for ");
    msg += get_material()->get_name() + " with ";
    msg += get_material()->get_structure();
    msg += " structure.";
    throw InitFailedException(msg);
  }
  add_submodel("semiconductor", semiconductor);

  //create kp hamiltonian for hole parameters
  kp_options["model"] = "6x6";

  bulk_ham = dynamic_cast<KPbulkHamiltonian*>(
      EFAbulkHamiltonian::create(get_material(), kp_options));

  if (bulk_ham == NULL)
    throw InitFailedException(string("Cannot create bulk hamiltonian for ")
        + get_material()->get_name());

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

  const Database& db = get_database();
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

void
DDsemiconductor::calculate_vb_bulk_states(const Tensor1& k_vector,
    std::vector<double>& eigenvalues, std::vector<libMesh::DenseVector<Complex> >& eigenvectors)
{

  eigenvalues.resize(6);
  eigenvectors.resize(6, libMesh::DenseVector<Complex>(6));

  bulk_ham->set_k_vector(k_vector);

  bulk_ham->calculate_Hamiltonian_k_par();

  if (strained) bulk_ham->apply_strain_and_potential(strain, 0.0);

  std::vector<std::vector<KPbulkHamiltonian::MatrixElement> >& Ham1 = bulk_ham->get_Hamiltonian();


  std::complex<double> ham6x6matrix[6*6];

  for (short i = 0; i < 6; i++)
    for (short j = i; j < 6; j++)
    {
      ham6x6matrix[i*6 + j] = 0.0;
      ham6x6matrix[j*6 + i] = Ham1[i][j].constant;
    }


  char jobs = 'V';
  char UPLO = 'U';
  int  N = 6;
  double eigvals[6];
  std::complex<double> WORK[11];
  int LWORK = 11;
  double RWORK[16];
  int info;

  zheev_(jobs, UPLO, N, ham6x6matrix, N, eigvals, WORK, LWORK, RWORK, info);
  if (info != 0 ) exit(1);

  for (short i = 0; i < 6; i++)
    eigenvalues[i] = eigvals[i] * Constants::Hartree;

  for (short i = 0; i < 6; i++)
    for (short j = 0; j < 6; j++)
      eigenvectors[i](j) = ham6x6matrix[i*6 + j];

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


    //cerr << "H = "  ;
    for (short i = 0; i < 6; i++)
      for (short j = 0; j < 6; j++)
	ham6x6matrix[i + j*6] = Ham1[i][j].constant;



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


void DDsemiconductor::calculate_inverse_mass(
    const vector<vector<KPbulkHamiltonian::MatrixElement> >& ham,
    const vector<libMesh::DenseVector<Complex> >::const_iterator first,
    const vector<libMesh::DenseVector<Complex> >::const_iterator last,
    map<ID, libMesh::RealTensor>& imasses)
{

  int N = distance(first, last);

  // calculate the matrix elements of the parameter matrices
  // with respect to the subspace basis
  // ordering here: xx yy zz xy xz yz
  vector<Complex*> Hkk(6);
  for (int n = 0; n < 6; ++n)
  {
    Hkk[n] = new Complex[N*N];
    for (int i = 0; i < N*N; i++)
      Hkk[n][i] = 0.0;
  }

  vector<libMesh::DenseVector<Complex> >::const_iterator it(first);
  for (int i = 0; it != last; ++i, ++it)
  {
    vector<libMesh::DenseVector<Complex> >::const_iterator jt(it);
    for (int j = i; jt != last; ++j, ++jt)
    {
      for (int a = 0; a < 6; a++)
      {
        for (int b = 0; b < 6; b++)
        {
          Complex val = conj((*it)(a)) * (*jt)(b);
          Hkk[0][j*N + i] +=  ham[a][b].quad[0][0] * val;
          Hkk[1][j*N + i] +=  ham[a][b].quad[1][1] * val;
          Hkk[2][j*N + i] +=  ham[a][b].quad[2][2] * val;
          Hkk[3][j*N + i] +=  (ham[a][b].quad[0][1] + ham[a][b].quad[1][0] +
              ham[a][b].quad[0][0] + ham[a][b].quad[1][1]) * val;
          Hkk[4][j*N + i] +=  (ham[a][b].quad[0][2] + ham[a][b].quad[2][0] +
              ham[a][b].quad[0][0] + ham[a][b].quad[2][2]) * val;
          Hkk[5][j*N + i] +=  (ham[a][b].quad[1][2] + ham[a][b].quad[2][1] +
              ham[a][b].quad[1][1] + ham[a][b].quad[2][2]) * val;
        }
      }
    }
  }


  Complex i(0,1);
  double is2 = 1/sqrt(2);
  double is3 = 1/sqrt(3);
  double is6 = 1/sqrt(6);
  double s23 = sqrt(2)/sqrt(3);
  Complex T[6][6];
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      T[i][j] = 0;
  T[0][0] = is2; T[0][1] = i*is2;
  T[1][3] = i*is2; T[1][4] = is2;
  T[2][0] = is6; T[2][1] = -i*is6; T[2][5] = s23;
  T[3][2] = -i*s23; T[3][3] = i*is6; T[3][4] = -is6;
  T[4][0] = -i*is3; T[4][1] = -is3; T[4][5] = i*is3;
  T[5][2] = is3; T[5][3] = is3; T[5][4] = i*is3;

  // solve the reduced eigensystem to obtain the inverse masses and the
  // appropriate expansion coefficients

  for (int k = 0; k < 6; ++k)
  {
    char jobs = 'V';
    char UPLO = 'U';
    double eigvals[6];
    std::complex<double> WORK[11];
    int LWORK = 11;
    double RWORK[16];
    int info;

    zheev_(jobs, UPLO, N, Hkk[k], N, eigvals, WORK, LWORK, RWORK, info);


    //
    // calculate the eigenvectors and identify band
    //

    std::vector<libMesh::DenseVector<Complex> > eigv;
    eigv.resize(N, libMesh::DenseVector<Complex>(6));

    // to avoid double counting of bands we keep track of already assigned ones
    set<int> identified_bands;
    for (int n = 0; n < N; n++)
    {
      // Hkk[k] contains the expansion coefficients with
      // respect to the subspace spanning eigenvectors
      for (int m = 0; m < N; m++)
        eigv[n].add(Hkk[k][n*N + m], *(first + m));

      // rotate according to calculation system
      // rotate according to direction of k
      rotate(eigv[n], static_cast<KDirection>(k));
      //rotate(eigv[n], get_material()->get_rotated_crystal());

      libMesh::DenseVector<Complex> vec(6);
      for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
          vec(i) += eigv[n](j) * conj(T[i][j]);

      vector<double> comp(6);
      int band = 0;
      for (int i = 0; i < 6; i++)
      {
        if (!identified_bands.count(i))
        {
          comp[i] = abs(vec(i));
          band = (comp[i] > comp[band]) ? i : band;
        }
      }
      identified_bands.insert(band);

      // insert inverse mass tensor component in map
      switch (k)
      {
        case X:
          imasses[band](0, 0) = 2.0 * eigvals[n];
          break;
        case Y:
          imasses[band](1, 1) = 2.0 * eigvals[n];
          break;
        case Z:
          imasses[band](2, 2) = 2.0 * eigvals[n];
          break;
        case XY:
          imasses[band](0, 1) = imasses[band](1, 0) = eigvals[n] -
            0.5 * (imasses[band](0, 0) + imasses[band](1, 1));
          break;
        case XZ:
          imasses[band](0, 2) = imasses[band](2, 0) = eigvals[n] -
            0.5 * (imasses[band](0, 0) + imasses[band](2, 2));
          break;
        case YZ:
          imasses[band](1, 2) = imasses[band](2, 1) = eigvals[n] -
            0.5 * (imasses[band](1, 1) + imasses[band](2, 2));
          break;

        default:
          throw RuntimeException("Unforeseen error: unknown k-direction.");

      }
    }

    delete Hkk[k];
  }


}


void DDsemiconductor::do_calculate_valence_band_extremum(void)
{


  vector<Tensor1> k_vector;
  k_vector.reserve(7);
  Tensor1 k;
  // Gamma
  k(1) = 0.0; k(2) = 0.0; k(3) = 0.0;
  k_vector.push_back(k);


  /*
  // solve the eigensystem in Gamma
  vector<double> eigval;
  vector<DenseVector<Complex> > eigvec;
  calculate_vb_bulk_states(k_vector[0], eigval, eigvec);


  // find groups of degenerate eigenvalues
  vector<unsigned int> groups(1, 1);
  for (size_t i = 1; i < eigval.size(); i++)
  {
    if (Utils::almost_equal::compare(eigval[i-1], eigval[i]))
      groups[groups.size() - 1]++;
    else
      groups.push_back(1);
  }

  // get the parameters for the k_i*k_j terms of the kp hamiltonian
  vector<vector<KPbulkHamiltonian::MatrixElement> > ham;
  bulk_ham->get_hamiltonian_without_k(ham);

  // the inverse masses, associated to bands
  map<ID, RealTensor> imasses;

  // iterate over the groups of degenerate states
  int end = 0;
  for (int g = 0; g < groups.size(); g++)
  {
    int N = groups[g];
    int first = end;
    end = first + N;

    cerr << endl << "E = " << eigval[first] << endl;
    calculate_inverse_mass(ham, eigvec.begin() + first, eigvec.begin() + end, imasses);

  }

  if (imasses.size() != 6)
    throw RuntimeException("Could not calculate all valence bands");

  {
  map<ID, RealTensor>::iterator it(imasses.begin());
  map<ID, RealTensor>::iterator end(imasses.end());
  for ( ; it != end; ++it)
  {
    cerr << "band = " << it->first << endl;
    cerr << it->second;
    double mass_DOS = std::pow(-1.0/it->second.det(),1.0/3.0);
    cerr << "mdos = " << mass_DOS << "  " << it->second.det() << endl << endl;
  }
  }
  */


  vector<DDsemiconductor::band_extremum> result;
  result.reserve(3);

  DDsemiconductor::band_extremum extremum;

  //-----------------------------------------------

  const Tensor2Gen& rotm = get_material()->get_rotated_crystal().RotMatrix;

  // [100]
  k(1) = k_max ; k(2) = 0; k(3) = 0;
  k_vector.push_back(rotm * k);

  // [010]
  k(1) = 0 ; k(2) = k_max; k(3) = 0;
  k_vector.push_back(rotm * k);

  // [001]
  k(1) = 0 ; k(2) = 0; k(3) = k_max;
  k_vector.push_back(rotm * k);

  // [110]
  k(1) = k_max; k(2) = k_max; k(3) = 0;
  k_vector.push_back(rotm * k);

  // [101]
  k(1) = k_max; k(2) = 0; k(3) = k_max;
  k_vector.push_back(rotm * k);

  // [011]
  k(1) = 0.0; k(2) = k_max; k(3) = k_max;
  k_vector.push_back(rotm * k);

  //--------------------------------------------------
  vector< vector<double> >  eigenvalue = calculate_vb_bulk_states(k_vector);

  Tensor2Sym imass;

  double Eh_kmax = Constants::Hartree * k_max * k_max;

  for (short ind = 0; ind < 3; ind++)
  {
    if (eigenvalue[0][ind*2] + energy_cutoff > eigenvalue[0][5])
    {
      extremum.degeneracy = 2;
      extremum.energy = eigenvalue[0][ind*2];

      imass(1,1) = (2.0 *(eigenvalue[0][ind*2] - eigenvalue[1][ind*2] )) / Eh_kmax;

      imass(2,2) = (2.0 *(eigenvalue[0][ind*2] - eigenvalue[2][ind*2] )) / Eh_kmax;

      imass(3,3) = (2.0 *(eigenvalue[0][ind*2] - eigenvalue[3][ind*2] )) / Eh_kmax;


      imass(2,1) = (eigenvalue[0][ind*2] - eigenvalue[4][ind*2]) / Eh_kmax
          - 0.5 * (imass(1,1) + imass(2,2));

      imass(3,1) = (eigenvalue[0][ind*2] - eigenvalue[5][ind*2]) / Eh_kmax
          - 0.5 * (imass(1,1) + imass(3,3));

      imass(3,2) = (eigenvalue[0][ind*2] - eigenvalue[6][ind*2]) / Eh_kmax
          - 0.5 * (imass(2,2) + imass(3,3));


      double imass_DOS = 0.0;
      double temp1, temp2;

      imass.invariants(&temp1, &temp2, &imass_DOS);

      if (imass_DOS <= 0.0)
      {
        Messages::warning("Negative valence band DOS mass, setting to 1.0");
	std::cout<<"band "<<ind<<std::endl;
	std::cout<<"Mass Tensor: "<<std::endl;
	std::cout<<imass(1,1)<<" "<<imass(2,1)<<" "<<imass(3,1)<<std::endl;
	std::cout<<imass(2,1)<<" "<<imass(2,2)<<" "<<imass(3,2)<<std::endl;
	std::cout<<imass(3,1)<<" "<<imass(3,2)<<" "<<imass(3,3)<<std::endl;

        imass_DOS = 1.0;
      }

      extremum.mass_DOS = std::pow(1.0/imass_DOS,1.0/3.0);


      result.push_back(extremum);
    }
  }

  valence_band = result;
}


