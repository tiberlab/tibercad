#include "UptWrapper.h"

//---------------------------------------------------------------------


UptWrapper::UptWrapper(){
  std::cout << "Constructing UPTIGHT instance... ";
    f77_upt_initsession(_handler);
    std::cout << "done." << std::endl;
    std::cout << "Received handler: ";
    for  (int ii = 0; ii < UPT_HSIZE; ++ii) {
      std::cout << _handler[ii] << " ";
    }
    std::cout << std::endl;
}


UptWrapper::~UptWrapper(){
  std::cout << "Destructing UPTIGHT instance... ";
    f77_upt_destructsession(_handler);
    std::cout << "done." << std::endl;
}


UptWrapper* UptWrapper::create()
{
  return new UptWrapper();
}



//!Assign simulation parameters to DFTB+ instance
void UptWrapper::fill_param(int verbose_lev, char *databasePath, char *workPath, 
                             char *gen_filename, char *gen_outname, int max_n_n, 
                             int harrison_flag, int relat_flag, int potential_flag, 
                             int optmat_flag, int poldir) {

  f77_upt_fillbasicparameters(_handler, verbose_lev, databasePath, workPath, gen_filename,
                                gen_outname, max_n_n, harrison_flag, relat_flag,
                                potential_flag, optmat_flag, poldir);


} 



//!Initialize UPT instance (allocations)
void UptWrapper::inituptight () {
    f77_upt_inituptight(_handler);
  }


void UptWrapper::add_potential(int nAtoms, double *potential)
{
	f77_upt_addpotential(_handler,nAtoms,potential);
}

//! add the k-points as a vector
void UptWrapper::add_kpoints(int numkp, double *k_vec)
{
	f77_upt_addkpoints(_handler,numkp,k_vec);
}



//! build ETB Hamiltonian with Uptight
void UptWrapper::compute_H () {
  f77_upt_createhamiltonian(_handler);
}



//! Lanczos diagonalization
void UptWrapper::lanczos_diag (int n_vb, int n_cb, double guess_vb, double guess_cb,
                                int min_iter, int long_iter, int max_iter, 
				double fast_tol, double long_tol, double ort_tol) {

  f77_upt_lanczosdiag (_handler, n_vb, n_cb, guess_vb, guess_cb, min_iter, long_iter,
		       max_iter, fast_tol, long_tol, ort_tol);


}







