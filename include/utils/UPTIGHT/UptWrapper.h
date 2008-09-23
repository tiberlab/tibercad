#ifndef _UPTWRAPPER_H_
#define _UPTWRAPPER_H_

#include "uptight.h"
#include <iostream>
#include <assert.h>
#include <complex>

//-----------------------------------------------------------------------

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


  //!Static method to create a Upt wrapper instance
  static UptWrapper* create(); 


  //!Function to fill Upt parameters:
  /*!
   * int  verbose_lev          : verbosity level
   * char(UPT_LC) databasePath : database path
   * char(UPT_LC) workPath     : work path
   * char(UPT_MC) gen_filename : gen filename to read structure
   * char(UPT_MC) gen_outname  : output name 
   * int max_n_n               : max n.n. order (only 1 works now)
   * int harrison_flag         : harrison scaling       (1 = true)
   * int relat_flag            : relativistic spin orbit
   * int potential_flag        : add external potential (1 = true)
   * int optmat_flag           : compute optical matrix (1 = true) 
   * int poldir                : polarized light x = 1, y = 2, z = 3 
   */
  void fill_param (int verbose_lev, char *databasePath, char *workPath, 
                   char *gen_filename, char *gen_outname, int max_n_n, 
                   int harrison_flag, int relat_flag, int potential_flag, 
		   int optmat_flag, int poldir);


  //!Initialize the Upt instance
  void inituptight ();


  //!Get UPTIGHT instance handler
  inline const int* get_handler(void){ return _handler; };


  //! Add an atom-projected potential to H
  //! Must be called after inituptight and before compute_H
  void add_potential(int nAtoms, double *potential);

	
  //! Add a list of k-points as a vector (3,numkp)
  void add_kpoints(int numkp, double *k_vec);

  
  //!Computes Hamiltonian (must be called after inituptight)
  void compute_H ();


  //!Function to perform Parravicini 2-step Lanczos
  /*!
   * int n_vb          : number of valence energy levels
   * int n_cb          : number of conduction energy levels
   * double guess_vb   : folding guess for valence
   * double guess_cb   : folding guess for conduction
   * int min_iter      : minimum number of iterations   (~2)
   * int long_iter     : number of long iterations      (~30)
   * int max_iter      : maximum number of iterations   (~3000)
   * double fast_tol   : tolerance for fast loop        (~1e-6)
   * double long_tol   : tolerance on long loop         (~1e-10)
   * double ort_tol    : orthogonality tolerance        (~1e-5)
   */
  void lanczos_diag (int n_vb, int n_cb, double guess_vb, double guess_cb,
                     int min_iter, int long_iter, int max_iter, 
                     double fast_tol, double long_tol, double ort_tol);

private:
  int _handler[UPT_HSIZE];

};


#endif
