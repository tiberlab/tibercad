#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include <vector>
#include <string>

//!this are the parameters to pass to SLEPc
struct SLEPCoptions
{
  unsigned int ev_number;

  std::string H_file_name;

  std::string S_file_name;

  std::string solver_type;

  double eps_tolerance;

  unsigned int eps_max_it;

  double st_ksp_rtol;

  std::string st_ksp_type;

  double spectrum_shift;
  
};


extern "C" 
{ 
  int  eig_value_problem_general(const SLEPCoptions& opt) ;
  void slepc_init();
  void slepc_done();
};




#endif
