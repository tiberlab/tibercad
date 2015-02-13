#ifndef _NEGFWRAPPER_H_
#define _NEGFWRAPPER_H_

#include "libnegf.h"
#include <iostream>
#include <assert.h>
#include <complex>
#include <vector>

//-----------------------------------------------------------------------
typedef std::complex<double> Complex; 

class NegfWrapper
{

public:

  //!Wrapper class for callings to libNEGF library (libnegf.so)

  //!Constructor
  /*!Assign an handler to libNEGF instance, transparent to programmer
   *
   */
  NegfWrapper();

  //! Destructor  
  ~NegfWrapper();

  //!Static method to create a Negf wrapper instance
  static NegfWrapper* create(); 

  //!Initialize the Negf instance
  void init();

  //!read H and S hamiltonians 
  int read_HS();

  //!read 'negf.in'
  int read_input();

  //!Compute current
  //! unitsOfH is a string like "H", "eV", "cm^-1"
  //! unitsOfJ is a string like "A"
  double current(std::string unitsOfH, std::string unitsOfJ);

  //!Compute charge density
  int density(std::vector<double>& density, std::string particle);

  //! Return the LDOS in matrix form
  void ldos(std::vector<double>& ldos, int esteps, int npoints);

  //!Set SC iteration
  void set_iteration(int iter);

  //!Set scratch folder path
  void set_scratch_path(std::string path);

  //!Set scratch folder path
  void set_output_path(std::string path);

  //!Clean the Negf instance variable space
  void clean_libnegf();

  //!Get UPTIGHT instance handler
  inline const int* get_handler(void){ return _handler; };

  //! Set verbosity level for the library screen output
  void set_verbose(int verbose_lev);

  //! used to set an integer value to label tunneling files
  void set_kpoint(int kpoint);

  //! used to set the reference contact 0= min mu; 1= max mu
  void set_reference(int minmax);

  //! Set what to do with outer parts 0=none; 1=upper block; 2=full
  void device_contact_dm(int outer);

  //! Get library version
  void get_version(void);

  //! used to set if print tunneling.dat 
  void set_write_tunn(bool flag);

  //! used to set if print LDOS.dat
  void set_write_ldos(bool flag);

  //! write partition info 
  void partition_info(void);
  
  //! set H via memcopy to libnegf container
  void set_H_csr(int nrow, std::vector<Complex >& A, std::vector<int>& JA,
                                                     std::vector<int>& IA );
  //! set S via memcopy to libnegf container
  void set_S_csr(int nrow, std::vector<Complex >& A, std::vector<int>& JA,
                                                     std::vector<int>& IA );
  //! Set S as and Identity matrix
  void set_S_id(int nrow);

  //! Print matrix
  void print_mat(void);

private:
  int _handler[NEGF_HSIZE];

};


#endif
