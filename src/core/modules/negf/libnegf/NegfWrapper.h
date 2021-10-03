#ifndef _NEGFWRAPPER_H_
#define _NEGFWRAPPER_H_

#include "libnegf.hpp"
#include "lnParams.h"

#include <complex>
#include <vector>

//-----------------------------------------------------------------------
typedef std::complex<double> Complex; 

//!Wrapper class for callings to libNEGF library (libnegf.so)
class NegfWrapper
{

public:

  typedef struct lnparams Parameters;

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

  //! Get the parameters
  void get_parameters(NegfWrapper::Parameters& params);

  //! Set the parameters
  void set_parameters(const NegfWrapper::Parameters& params);

  //!Compute current
  //! unitsOfH is a string like "H", "eV", "cm^-1"
  //! unitsOfJ is a string like "A"
  double current(std::string unitsOfH, std::string unitsOfJ);

  //!Compute charge density
  void density(std::vector<double>& density, std::string particle);

  //! Return the LDOS in matrix form
  void ldos(std::vector<double>& ldos, int& esteps, int& npoints);

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

  //! Initialized the contacts
  void init_contacts(int n_cont);

  //! Initialize structure
  void init_structure(int ncont, int* surfstart, 
    int* surfend, int* contend, int npl, int* plend, int* cblk);

  //! Get library version
  void get_version(void);

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
