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
			    char *outPath,
			    char *gen_filename, char *gen_outname, char *sparse_fmt, 
			    int max_n_n, bool harrison, bool relat, bool potential,
			    bool optmat, int poldir, double *c_axis, bool check_bondmap) {

  int harrison_flag, relat_flag, potential_flag, optmat_flag, check_nn;
  
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

  f77_upt_fillbasicparameters(_handler, verbose_lev, databasePath, workPath, outPath,
			      gen_filename, gen_outname, sparse_fmt, max_n_n, 
			      harrison_flag, relat_flag, potential_flag, optmat_flag, 
			      poldir, c_axis[0], c_axis[1], c_axis[2], check_nn);

}



//!Initialize UPT instance (allocations)
void UptWrapper::inituptight () {
    f77_upt_inituptight(_handler);
  }

//!destroy container variables (allocations)
void UptWrapper::cleanuptight() 
{
        f77_upt_destructuptight(_handler);
}


void UptWrapper::add_potential(std::vector<double>& potential)
{
  
  f77_upt_addpotential(_handler,potential.size(),&potential.front());
}

//! add the k-points as a vector
void UptWrapper::add_kpoints(int numkp, double *k_vec)
{
	f77_upt_addkpoints(_handler,numkp,k_vec);
}



//! build ETB Hamiltonian with Uptight
void UptWrapper::compute_H(){
  f77_upt_createhamiltonian(_handler);
}

//! build ETB Hamiltonian with Uptight
void UptWrapper::compute_P_matrix(int poldir){

  f77_upt_setpmatrix(_handler, 1, poldir);
  
  f77_upt_createhamiltonian(_handler);

  f77_upt_setpmatrix(_handler, 0, poldir);

}



//! Lanczos diagonalization
void UptWrapper::lanczos_diag(int n_vb, int n_cb, double guess_vb, double guess_cb,
                                int min_iter, int long_iter, int max_iter,
				double fast_tol, double long_tol, double ort_tol) {

  f77_upt_lanczosdiag(_handler, n_vb, n_cb, guess_vb, guess_cb, min_iter, long_iter,
		       max_iter, fast_tol, long_tol, ort_tol);


}

void UptWrapper::set_num_states(int n_vb, int n_cb)
{
  f77_upt_set_num_states(_handler,n_vb,n_cb);
}


//! get ETB Hamiltonian size
int UptWrapper::get_H_dim(void) {
  int hdim;	
  f77_upt_get_hamildim(_handler,hdim);
  return hdim;
}

//! write eigenstates on file
void UptWrapper::write_states() {	
  f77_upt_write_states(_handler);
}

//! read eigenstates from file
void UptWrapper::read_old_states(void) {	
  f77_upt_read_states(_handler);
}


//! get computed states 
void UptWrapper::get_states(int num_ev, int hdim,
			     double* eigenvals, double *eigenvec_re, double *eigenvec_im) {

  f77_upt_get_states (_handler, num_ev, hdim, eigenvals, eigenvec_re, eigenvec_im);


}


std::complex<double> UptWrapper::get_matel(int i, int j)
{

  std::complex<double> matel;
  double matel_re, matel_im;

  f77_upt_get_matel(_handler,i,j,matel_re,matel_im);

  matel = std::complex<double>(matel_re,matel_im);

  //std::cerr<<matel<<std::endl;

  return matel;
}

void UptWrapper::get_ion_numorbitals(std::vector<int>& ion_block_vector)
{
  //int nAtoms = ion_block_vector.size(); 
  //int* p_ion_block_vec = new int[nAtoms];

  f77_upt_get_ion_numorbitals(_handler,&ion_block_vector.front());

  //for (int j = 0; j < nAtoms; j++) 
  //  ion_block_vector[j] = p_ion_block_vec[j];

  //delete [] p_ion_block_vec;

}


void UptWrapper::set_verbose(int verbose_lev)
{

  f77_upt_set_verbosity(_handler,verbose_lev);

}



