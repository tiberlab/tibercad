#include "UptWrapper.h"


//---------------------------------------------------------------------


UptWrapper::UptWrapper(){
    std::cout << "\nCreating UPTIGHT instance... ";
    upt_initsession_(_handler);
    std::cout << _handler << std::endl;

    std::cout << std::endl;
    upt_getversion_(_handler);
}


UptWrapper::~UptWrapper(){
    upt_destructuptight_(_handler);
    upt_destructsession_(_handler);
}


UptWrapper* UptWrapper::create()
{
  return new UptWrapper();
}


void UptWrapper::set_paths(const char* databasePath, 
                          const char* workPath, const char* outPath) 
{

   upt_set_paths_(_handler, databasePath, workPath, outPath);
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

  upt_fillbasicparameters_(_handler, verbose_lev, 
			      gen_filename, gen_outname, sparse_fmt, max_n_n, 
			      harrison_flag, relat_flag, potential_flag, optmat_flag, 
			      poldir, c_axis[0], c_axis[1], c_axis[2], check_nn, 
			      dg_scale, dg_onsite, hybrid );

}

void UptWrapper::set_solver_flag(int flag)
{ 
  upt_solver_flag_(_handler, flag);
}

void UptWrapper::set_output(int format, double step)
{
  upt_setoutput_(_handler, format, step);
}

int UptWrapper::inituptight() 
{
    int hsize;
    upt_gethandlersize_(hsize);

    if(hsize!=UPT_HSIZE){
      return 1;      
    }
    else
    {
      upt_inituptight_(_handler);
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
      upt_setmpicomm_(_handler, MPI_Comm_c2f(comm));
      return 0;
    }
}

void UptWrapper::get_version()
{
    upt_getversion_(_handler);
}


//!destroy container variables (allocations)
void UptWrapper::cleanuptight() 
{
    upt_destructuptight_(_handler);
}


void UptWrapper::add_potential(std::vector<double>& potential)
{
int size = potential.size();  
  upt_addpotential_(_handler,size,&potential.front());
}


void UptWrapper::clear_potential(void)
{
  
  upt_erasepotential_(_handler);
}



void UptWrapper::set_kpoint(double *k_vec)
{
  upt_setkpoint_(_handler, k_vec);
}



void UptWrapper::compute_H(char* sprs_fmt){	
  upt_createhamiltonian_(_handler, sprs_fmt);
}

void UptWrapper::print_H(void){	
  upt_printhamiltonian_(_handler);
}


void UptWrapper::compute_P_matrix(int poldir, char* sprs_fmt)
{
  upt_setpmatrix_(_handler, 1, poldir);
  upt_createhamiltonian_(_handler, sprs_fmt);
  upt_setpmatrix_(_handler, 0, poldir);
}


void UptWrapper::lanczos_diag(int st_vb, int st_cb, int n_vb, int n_cb, double guess_vb, double guess_cb,
                                int min_iter, int long_iter, int max_iter,
				                        double fast_tol, double long_tol, double ort_tol,
				                      int dynamic, double bitoff) {

  upt_lanczosdiag_(_handler, st_vb, st_cb, n_vb, n_cb, guess_vb, guess_cb, min_iter, long_iter,
		       max_iter, fast_tol, long_tol, ort_tol, dynamic, bitoff);


}


void UptWrapper::jacobidavidson(int st_cb, int st_vb, int n_vb, int n_cb, double guess_vb, double guess_cb,
      double long_tol)
{
  upt_jd_diag_(_handler, st_cb, st_vb, n_vb, n_cb, guess_vb, guess_cb, long_tol);
}


void UptWrapper::feast(double emin, double emax, int m0) {

  upt_feastsolver_(_handler, emin, emax, m0);

}

void UptWrapper::lapack(int n_vb, int n_cb, double guess_vb, double guess_cb) {

  upt_lapacksolver_(_handler, n_vb, n_cb, guess_vb, guess_cb);

}

void UptWrapper::set_num_states(int n_vb, int n_cb)
{
  upt_set_num_states_(_handler,n_vb,n_cb);
}


int UptWrapper::get_H_dim(void) {
  int hdim;	
  upt_get_hamildim_(_handler,hdim);
  return hdim;
}

int UptWrapper::get_H_nnz(void) {
  int hdim;	
  upt_get_hamilnnz_(_handler,hdim);
  return hdim;
}

int UptWrapper::get_H_row_size(int row) {
   int size;
   upt_get_hamil_rowsize_(_handler,row,size);
   return size;
}

void UptWrapper::get_H_row(int row, int* colind, Complex* vals) {
   upt_get_hamil_row_(_handler,row,colind,vals);
}

void UptWrapper::set_statefile(const char* filename)
{
  upt_setstatefile_(_handler, filename);
}

void UptWrapper::write_states() {
  upt_write_states_(_handler);
}

void UptWrapper::set_outpath(const char* outpath) {
  upt_setoutpath_(_handler, outpath); 
}

void UptWrapper::set_workpath(const char* workpath) {
  upt_setworkpath_(_handler, workpath); 
}

void UptWrapper::read_old_states(char* load_path, int& nev, int& nec) {	
  //upt_setloadpath(_handler, load_path);
  upt_read_states_(_handler, load_path, nev, nec);
}


void UptWrapper::get_states(int num_ev, int hdim, double* eigenvals, 
                            Complex* eigenvec, int* particles) {

  upt_get_states_ (_handler, num_ev, hdim, eigenvals, eigenvec, particles);


}


void UptWrapper::set_state(int num_ev, int id, int hdim, double eigval,
    const std::vector<Complex>& state, int particle)
{
  upt_set_state_(_handler, num_ev, id+1, hdim, eigval, state.data(), particle);
}


Complex UptWrapper::get_matel(int i, int j)
{

  Complex matel;

  upt_get_matel_(_handler,i,j,matel);

  return matel;
}

void UptWrapper::get_ion_numorbitals(std::vector<int>& ion_block_vector)
{
  upt_get_ion_numorbitals_(_handler, &ion_block_vector.front());
}


void UptWrapper::get_ion_orbitals(int i, std::vector<int>& orbitals)
{
  upt_get_ion_orbitals_(_handler, i, &orbitals.front());
}


void UptWrapper::set_verbose(int verbose_lev)
{

  upt_set_verbosity_(_handler,verbose_lev);

}


void UptWrapper::set_strain(std::vector<double>& e_xx, std::vector<double>& e_yy, 
                            std::vector<double>& e_zz)
{
  int size = e_xx.size();	
  upt_setstrain_(_handler, size, &e_xx.front(), &e_yy.front(), &e_zz.front());
}
 

void UptWrapper::get_H_csr(int nrow, char fmt, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{

  upt_get_csr_hamiltonian_(_handler,nrow,fmt,&A.front(),&JA.front(),&IA.front());

}

void UptWrapper::set_H_csr(int nrow, char fmt, std::vector<Complex >& A, 
                           std::vector<int>& JA, std::vector<int>& IA)
{

  upt_set_csr_hamiltonian_(_handler,nrow,fmt,&A.front(),&JA.front(),&IA.front());

}

void UptWrapper::complex_test(double& re, double& im, Complex& zz)
{
  complex_test_(re,im,zz);
}

double UptWrapper::real_test()
{
  double re;
  real_test_(re);
  return re;
}
