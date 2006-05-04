#ifndef _SBBULKHAMILTONIAN_H_
#define _SBBULKHAMILTONIAN_H_
//! A class that builds single band Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"

class SBbulkHamiltonian : public EFAbulkHamiltonian
{

 public:
  //! default constructor
  SBbulkHamiltonian(void);

  //! constructor
  /*!
    \param band_edge band edge energy [eV]
    \param imass 1/m tensor [m0]
  */
  SBbulkHamiltonian(double band_edge, Tensor2Sym& imass);


  //! set edge energy
  /*!
    \param energy band edge energy [eV]
  */
  void set_band_edge_energy(double energy);


  void set_diag_mass_tensor(double m_xx, double m_yy, double m_zz);

  virtual void calculate_Hamiltonian_k_par(void);

 

  virtual void calculate_Hamiltonian_gen(void); 

  //! Applies ONLY potential. Strain in applied only in derived classes
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

 private:

  double edge;

  Tensor2Sym imass;

  MatrixElement single_band_ham;

};

#endif
