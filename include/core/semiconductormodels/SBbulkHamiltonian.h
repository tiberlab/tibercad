#ifndef _SBBULKHAMILTONIAN_H_
#define _SBBULKHAMILTONIAN_H_
//! A class that builds single band Hamiltonian

#include <complex>
#include <vector>
#include "tensor.h"
#include "EFAbulkHamiltonian.h"

class DDsemiconductor;


class SBbulkHamiltonian : public EFAbulkHamiltonian
{

 public:
  //! default constructor
  SBbulkHamiltonian(void);

  //! destructor
  ~SBbulkHamiltonian(void);

  //! constructor
  /*!
    \param band_edge band edge energy [eV]
    \param imass 1/m tensor [m0]
  */
  SBbulkHamiltonian(double band_edge, Tensor2Sym& imass);


  //!copy constructor
  SBbulkHamiltonian(const SBbulkHamiltonian& model);


  //! set edge energy
  /*!
    \param energy band edge energy [eV]
  */
  void set_band_edge_energy(double energy);


  void set_diag_mass_tensor(double m_xx, double m_yy, double m_zz);

  void calculate_Hamiltonian_k_par(void);

 

  void calculate_Hamiltonian_gen(void); 


  //! Applies ONLY potential. Strain in applied only in derived classes
  virtual void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

   
  void read_database(const Dummy&);

  //! \deprecated { Create parameters for an alloy }
  /*!
   * \deprecated { This method will live as long as the database is
   * not used yet.}
   */

  void  build_alloy(const std::string& component2,
			   const std::string& bowing_params, double content);

  void set_data_file(const std::string& filename)
      { _filename = filename; };

  
  //! here does nothing but may be  overloaded in the derived classes.
  virtual void  calculate_edge_and_mass() = 0;
  



 protected:

  double edge;

  Tensor2Sym imass;

  MatrixElement single_band_ham;

  std::string _filename;

 
  //! a pointer to a semiconductor that contains parameters
  DDsemiconductor* semiconductor;
  
 private:

 
};

#endif
