#ifndef _EFAbulkHamiltonian_H_
#define _EFAbulkHamiltonian_H_
//! A general class for envelope function bulk Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"

class EFAbulkHamiltonian
{

 public:

  typedef std::complex<double> Complex;

  struct MatrixElement
  {
    std::complex<double> constant;
    std::complex<double> linear_left[3];
    std::complex<double> linear_right[3];
    std::complex<double> quad[3][3];
    
  };

  EFAbulkHamiltonian();

  void set_k_vector(const double kvector[3]);
 
  void set_k_vector(Tensor1 k_in); 

  void set_rotation_matrix(Tensor2Gen& rotmatrix);

  //! calculate model Hamiltonian without application of k||
  virtual void calculate_Hamiltonian_gen(void) = 0; 
 
  //! apply k|| to the Hamiltonian
  virtual void calculate_Hamiltonian_k_par(void) = 0;

  std::vector< std::vector<MatrixElement > >& get_Hamiltonian(void); 

  //!destructor
  virtual ~EFAbulkHamiltonian() {};
  
 protected:
  

  //!Hamiltonian in k representation
  std::vector< std::vector<MatrixElement > > Hamiltonian;

  //! k-vector in simualtion system
  double k_vector[3];


  //! Hartree energy in eV
  static const double Hartree = 27.2113961;


  //!rotation matrix
  double  rot_matrix[3][3];
  
  //!rotates tensor of rank 1 from crystal to calculation system
  void    rotate_linear(std::complex<double> *vector);

  //!rotates tensor of rank 2 from crystal to calculation system
  void    rotate_quad(std::complex<double> matrix[][3]);



 private:
  
  MatrixElement Ham;
};

#endif
