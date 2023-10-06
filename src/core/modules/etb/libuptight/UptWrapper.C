#include "UptWrapper.h"


//---------------------------------------------------------------------


UptWrapper::UptWrapper(){
    std::cout << "\nCreating UPTIGHT instance... ";
    f77_upt_initsession(_handler);
    std::cout << _handler << std::endl;

    std::cout << std::endl;
    f77_upt_getversion(_handler);
}


UptWrapper::~UptWrapper(){
    f77_upt_destructuptight(_handler);
    f77_upt_destructsession(_handler);
}


UptWrapper* UptWrapper::create()
{
  return new UptWrapper();
}


void UptWrapper::set_paths(const char* databasePath, 
                          const char* workPath, const char* outPath) 
{

   f77_upt_set_paths(_handler, databasePath, workPath, outPath);
}

void UptWrapper::fill_param(int verbose_lev,  
			    char *gen_filename, char *gen_outname, char *sparse_fmt, 
			    int max_n_n, bool harrison, bool relat, bool potential,
			    bool optmat, int poldir, double *c_axis, bool check_bondmap,
			    double dg_scale, double dg_onsite, bool hybrid_passivation)
{

  int harrison_flag, relat_flag, potential_flag, optmat_flag, check_nn, hybrid;
  
  harrison_flag=0;
  if(harrison){harrison_flag=1;}
  relat_flag=0;
  if(relat){relat_flag=1;}
  potential_flag=0;
  if(potential){potential_flag=1;}
  optmat_flag=0;
  if(optmat){optmat_flag=1;}
  check_nn=0;
  if(check_bondmap){check_nn=1;}
  hybrid=0;
  if (hybrid_passivation){hybrid=1;}

  f77_upt_fillbasicparameters(_handler, verbose_lev, 
			      gen_filename, gen_outname, sparse_fmt, max_n_n, 
			      harrison_flag, relat_flag, potential_flag, optmat_flag, 
			      poldir, c_axis[0], c_axis[1], c_axis[2], check_nn, 
			      dg_scale, dg_onsite, hybrid );

}

void UptWrapper::set_solver_flag(int flag)
{ 
  f77_upt_solver_flag(_handler, flag);
}

void UptWrapper::set_output(int format, double step)
{
  f77_upt_setoutput(_handler, format, step);
}

int UptWrapper::inituptight() 
{
    int hsize;
    f77_upt_gethandlersize(hsize);

    if(hsize!=UPT_HSIZE){
      return 1;      
    }
    else
    {
      f77_upt_inituptight(_handler);
      return 0;
    }
}

int UptWrapper::set_mpi_comm(MPI_Comm comm)
{
    // 4/5/2016 simply casting gives an error, we need to correctly
    // convert the handle

    // checks that MPI_comm can be casted to integer
    //int icomm = (int) comm;
    //if (icomm != comm) 
    //{
    //  return 1;
    //}
    //else
    {
      f77_upt_setmpicomm(_handler, MPI_Comm_c2f(comm));
      return 0;
    }
}

void UptWrapper::get_version()
{
    f77_upt_getversion(_handler);
}


//!destroy container variables (allocations)
void UptWrapper::cleanuptight() 
{
    f77_upt_destructuptight(_handler);
}


void UptWrapper::add_potential(std::vector<double>& potential)
{
int size = potential.size();  
  f77_upt_addpotential(_handler,size,&potential.front());
}


void UptWrapper::clear_potential(void)
{
  
  f77_upt_erasepotential(_handler);
}



void UptWrapper::set_kpoint(double *k_vec)
{
  f77_upt_setkpoint(_handler, k_vec);
}



void UptWrapper::compute_H(char* sprs_fmt){	
  f77_upt_createhamiltonian(_handler, sprs_fmt);
}

void UptWrapper::print_H(void){	
  f77_upt_printhamiltonian(_handler);
}


void UptWrapper::compute_P_matrix(int poldir, char* sprs_fmt)
{
  f77_upt_setpmatrix(_handler, 1, poldir);
  f77_upt_createhamiltonian(_handler, sprs_fmt);
  f77_upt_setpmatrix(_handler, 0, poldir);
}


void UptWrapper::lanczos_diag(int st_vb, int st_cb, int n_vb, int n_cb, double guess_vb, double guess_cb,
                                int min_iter, int long_iter, int max_iter,
				                        double fast_tol, double long_tol, double ort_tol,
				                      int dynamic, double bitoff) {

  f77_upt_lanczosdiag(_handler, st_vb, st_cb, n_vb, n_cb, guess_vb, guess_cb, min_iter, long_iter,
		       max_iter, fast_tol, long_tol, ort_tol, dynamic, bitoff);


}


void UptWrapper::jacobidavidson(int st_cb, int st_vb, int n_vb, int n_cb, double guess_vb, double guess_cb,
      double long_tol)
{
  f77_upt_jd_diag(_handler, st_cb, st_vb, n_vb, n_cb, guess_vb, guess_cb, long_tol);
}


void UptWrapper::feast(double emin, double emax, int m0) {

  f77_upt_feastsolver(_handler, emin, emax, m0);

}

void UptWrapper::lapack(int n_vb, int n_cb, double guess_vb, double guess_cb) {

  f77_upt_lapacksolver(_handler, n_vb, n_cb, guess_vb, guess_cb);

}

void UptWrapper::set_num_states(int n_vb, int n_cb)
{
  f77_upt_set_num_states(_handler,n_vb,n_cb);
}


int UptWrapper::get_H_dim(void) {
  int hdim;	
  f77_upt_get_hamildim(_handler,hdim);
  return hdim;
}

int UptWrapper::get_H_nnz(void) {
  int hdim;	
  f77_upt_get_hamilnnz(_handler,hdim);
  return hdim;
}

int UptWrapper::get_H_row_size(int row) {
   int size;
   f77_upt_get_hamil_rowsize(_handler,row,size);
   return size;
}

void UptWrapper::get_H_row(int row, int* colind, Complex* vals) {
   f77_upt_get_hamil_row(_handler,row,colind,vals);
}

void UptWrapper::set_statefile(const char* filename)
{
  f77_upt_setstatefile(_handler, filename);
}

void UptWrapper::write_states() {
  f77_upt_write_states(_handler);
}

void UptWrapper::set_outpath(const char* outpath) {
  f77_upt_setoutpath(_handler, outpath); 
}

void UptWrapper::set_workpath(const char* workpath) {
  f77_upt_setworkpath(_handler, workpath); 
}

void UptWrapper::read_old_states(char* load_path, int& nev, int& nec) {	
  //f77_upt_setloadpath(_handler, load_path);
  f77_upt_read_states(_handler, load_path, nev, nec);
}


void UptWrapper::get_states(int num_ev, int hdim, double* eigenvals, 
                            Complex* eigenvec, int* particles) {

  f77_upt_get_states (_handler, num_ev, hdim, eigenvals, eigenvec, particles);


}


void UptWrapper::set_state(int num_ev, int id, int hdim, double eigval,
    const std::vector<Complex>& state, int particle)
{
  f77_upt_set_state(_handler, num_ev, id+1, hdim, eigval, state.data(), particle);
}


Complex UptWrapper::get_matel(int i, int j)
{

  Complex matel;
  //double matel_re, matel_im;

  f77_upt_get_matel(_handler,i,j,matel);

  //matel = Complex(matel_re,matel_im);

  //std::cerr<<matel<<std::endl;

  return matel;
}

void UptWrapper::get_ion_numorbitals(std::vector<int>& ion_block_vector)
{
  //int nAtoms = ion_block_vector.size(); 
  //int* p_ion_block_vec = new int[nAtoms];

  f77_upt_get_ion_numorbitals(_handler, &ion_block_vector.front());

  //for (int j = 0; j < nAtoms; j++) 
  //  ion_block_vector[j] = p_ion_block_vec[j];

  //delete [] p_ion_block_vec;

}


void UptWrapper::get_ion_orbitals(int i, std::vector<int>& orbitals)
{
  f77_upt_get_ion_orbitals(_handler, i, &orbitals.front());
}


void UptWrapper::set_verbose(int verbose_lev)
{

  f77_upt_set_verbosity(_handler,verbose_lev);

}


void UptWrapper::set_strain(std::vector<double>& e_xx, std::vector<double>& e_yy, 
                            std::vector<double>& e_zz)
{
  int size = e_xx.size();	
  f77_upt_setstrain(_handler, size, &e_xx.front(), &e_yy.front(), &e_zz.front());
}
 

void UptWrapper::get_H_csr(int nrow, char fmt, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{

  f77_upt_get_csr_hamiltonian(_handler,nrow,fmt,&A.front(),&JA.front(),&IA.front());

}

void UptWrapper::set_H_csr(int nrow, char fmt, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{

  f77_upt_set_csr_hamiltonian(_handler,nrow,fmt,&A.front(),&JA.front(),&IA.front());

}

void UptWrapper::complex_test(double& re, double& im, Complex& zz)
{
  f77_complex_test(re,im,zz);
}

double UptWrapper::real_test()
{
  double re;
  f77_real_test(re);
  return re;
}
