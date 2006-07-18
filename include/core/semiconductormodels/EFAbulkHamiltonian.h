#ifndef _EFAbulkHamiltonian_H_
#define _EFAbulkHamiltonian_H_
//! A general class for envelope function bulk Hamiltonian

#include <complex>
#include <vector>
#include <map>
#include "tensor.h"
#include "PhysicalProperties.h"
class EFAbulkHamiltonian: public PhysicalProperties
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

  //! apply strain and potential to the EFA Hamiltonian
  /*!
    \param strain_crystal strain tensor in crystal system
    \param el_potentia electric potential [V]
  */
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential) = 0;

  std::vector< std::vector<MatrixElement > >& get_Hamiltonian(void); 

  //!destructor
  virtual ~EFAbulkHamiltonian() {};
  //----------------------------------------------------------------------------------
  //virtual methods of  PhysicalProperties
  virtual void read_database (const Dummy &db);

  virtual void read_database_bowing_parameters (const Dummy &db);
  
  virtual void set_properties_alloy (const PhysicalProperties *prop_comp1, const PhysicalProperties *prop_comp2, double molar_fraction);
  //---------------------------------------------------------------------------------

  const std::map<short, short>&  get_kp_bands_map(void) const; 

 protected:
  
  


  //!result Hamiltonian in k representation 
  std::vector< std::vector<MatrixElement > > Hamiltonian;


  //! Hamiltonian in k representation that is used by apply_strain_and_potential
  std::vector< std::vector<MatrixElement > > Hamiltonian_without_strain_pot;

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

  //!numbers of the bands in a 8x8 k.p basis
  /*
    0,1 - conduction bands;
    2-7 - valence bands
  */
  std::vector<short> kp_bands;


  //!map between band numbers  
  std::map <short, short> kp_bands_map;

 private:
  
  // MatrixElement Ham;
};

#endif
