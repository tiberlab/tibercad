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
  int init_libnegf();

  //!Clean the Negf instance variable space
  void clean_libnegf();

  //!Get UPTIGHT instance handler
  inline const int* get_handler(void){ return _handler; };

  //! Set verbosity level for the library screen output
  void set_verbose(int verbose_lev);

  //! Get library version
  void get_version(void);

private:
  int _handler[NEGF_HSIZE];

};


#endif
