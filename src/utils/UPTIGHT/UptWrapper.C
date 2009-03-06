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



//!Assign simulation parameters to UPT instance
void UptWrapper::fill_param(int verbose_lev, char *databasePath, char *workPath,
                             char *gen_filename, char *gen_outname, int max_n_n,
                             bool harrison, bool relat, bool potential,
                             bool optmat, int poldir, double *c_axis) {

  int harrison_flag, relat_flag, potential_flag, optmat_flag;
  
  harrison_flag=0;
  if(harrison){harrison_flag=1;}
  relat_flag=0;
  if(relat){relat_flag=1;}
  potential_flag=0;
  if(potential){potential_flag=1;}
  optmat_flag=0;
  if(optmat){optmat_flag=1;}

  f77_upt_fillbasicparameters(_handler, verbose_lev, databasePath, workPath, gen_filename,
                                gen_outname, max_n_n, harrison_flag, relat_flag,
                                potential_flag, optmat_flag, poldir, c_axis[0], c_axis[1],
			        c_axis[2]);

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







