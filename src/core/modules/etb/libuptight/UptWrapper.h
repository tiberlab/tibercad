#ifndef _UPTWRAPPER_H_
#define _UPTWRAPPER_H_

#include "uptight.h"
#include <iostream>
#include <assert.h>
#include <complex>
#include <vector>

#include <mpi.h>
//-----------------------------------------------------------------------
typedef std::complex<double> Complex; 

class UptWrapper
{

public:

  //!Wrapper class for callings to DFTB+ library (libdftbp.so)



  //!Constructor
  /*!Assign an handler to DFTB+ instance, transparent to programmer
   *
   */
  UptWrapper();


  //! Destructor  
  ~UptWrapper();


  //! Static method to create a Upt wrapper instance
  static UptWrapper* create(); 

  //! Initialize the MPI communicator for the library
  int set_mpi_comm(MPI_Comm comm);

  //!Function to fill Upt parameters:
  /*!
   * int  verbose_lev          : verbosity level
   * char(UPT_LC) databasePath : database path
   * char(UPT_LC) workPath     : work path
   * char(UPT_MC) gen_filename : gen filename to read structure
   * char(UPT_MC) gen_outname  : output name 
   * int max_n_n               : max n.n. order (only 1 works now)
   * bool harrison_flag        : harrison scaling      
   * bool relat_flag           : relativistic spin orbit
   * bool potential_flag       : add external potential
   * bool optmat_flag          : compute optical matrix
   * int poldir                : polarized light x = 1, y = 2, z = 3 
   * c_axis                    : set direction of c-axis in wurtzites
   * check_bondmap             : do check bondmap consistency
   * dg_scale                  : scaling factor for dg-bonds
   * dg_onsite                 : scaling factor for dg-onsite energies
   * hybrid_passivation        : perform Boykin-type passivation
   */
  void fill_param(int verbose_lev, char *gen_filename, char *gen_outname, char *sparse_fmt, 
		   int max_n_n, bool harrison_flag, bool relat_flag, 
		   bool potential_flag, bool optmat_flag, int poldir, 
		   double *c_axis, bool check_bondmap, 
                   double dg_scale, double dg_onsite, bool hybrid_passivation);

  //! Set parameters for eigenstates output format 
  //! out_format == 1 => JVXL
  //! out_format == 2 => CUBE
  void set_output(int out_format, double scale);

  void set_paths(const char* defaultPath, const char* databasePath, 
                 const char* workPath, const char* outPath); 

  void set_workpath(const char* workPath);
  
  void set_outpath(const char* outPath);

  void set_statefile(const char* filename);

  //! Set verbosity level for the library screen output
  void set_verbose(int verbose_lev);

  //!Initialize the Upt instance 
  //!Read geometry and build n.n. table (call after setting globals and fill_param) 
  //! Get library version
  void get_version(void);

  //!Initialize the Upt instance
  int inituptight();

  //!Clean the Upt instance variable space
  void cleanuptight();

  //!Get UPTIGHT instance handler
  inline const int* get_handler(void){ return _handler; };

  //! Add an atom-projected potential to H
  //! Must be called after inituptight and before compute_H
  void add_potential(std::vector<double>& potential);


  //! set atom-projected strain
  void set_strain(std::vector<double>& e_xx, std::vector<double>& e_yy, std::vector<double>& e_zz);


  void clear_potential(void);
	
  //! Set a k-point, k_vec(3)
  void set_kpoint(double *k_vec);

  
  //!Computes Hamiltonian (must be called after inituptight)
  void compute_H(char* sprs_fmt);

  //!Computes Hamiltonian (must be called after inituptight)
  void print_H(void);

  //!Computes P-matrix for optics (must be called after inituptight)
  void compute_P_matrix(int poldir, char* sprs_fmt);


  //!Function to perform Parravicini 2-step Lanczos
  /*!
   * int n_vb          : number of valence energy levels
   * int n_cb          : number of conduction energy levels
   * double guess_vb   : folding guess for valence
   * double guess_cb   : folding guess for conduction
   * int min_iter      : minimum number of iterations   (~2)
   * int long_iter     : number of long iterations      (~30)
   * int max_iter      : maximum number of iterations   (~10000)
   * double fast_tol   : tolerance for fast loop        (~1e-1)
   * double long_tol   : tolerance on long loop         (~1e-10)
   * double ort_tol    : orthogonality tolerance        (~1e-5)
   * int dynamic       : new way using dynamical shift
   */
  void lanczos_diag (int st_vb, int st_cb, int n_vb, int n_cb, double guess_vb, double guess_cb,
                     int min_iter, int long_iter, int max_iter, 
                     double fast_tol, double long_tol, double ort_tol,
		                 int dynamic);

  //! Interface for Jacobi-Davidson
  void jacobidavidson(int st_cb, int st_vb, int n_vb, int n_cb, double guess_vb, double guess_cb,
      double long_tol);


  // Solver flag. currently 0: CPU solver; 1: GPU; 2: GPU-split version 
  void set_solver_flag(int flag);

  void set_num_states(int n_vb, int n_cb);

  int get_H_dim(void);

  int get_H_nnz(void);

  int get_H_row_size(int row);

  void get_H_row(int row, int* colind, Complex* vals);
  
  void get_H_csr(int nrow, char fmt, std::vector<Complex >& A, std::vector<int>& JA,
                                                     std::vector<int>& IA );
  
  void set_H_csr(int nrow, char fmt, std::vector<Complex >& A, std::vector<int>& JA,
                                                     std::vector<int>& IA );

  void write_states(void);

  void read_old_states(char* path, int& nev, int& nec);

  void get_states(int num_ev, int hdim, double* eigenvals,
                  Complex* states, int* particles); 

  Complex get_matel(int i, int j);

  void get_ion_numorbitals(std::vector<int>& ion_block_vector);


  void complex_test(double& re, double& im, Complex& zz);

  double real_test();

  void feast(double emin, double emax, int m0) ;

private:
  int _handler[UPT_HSIZE];

};


#endif
