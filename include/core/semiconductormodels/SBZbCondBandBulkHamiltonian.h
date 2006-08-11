#ifndef _SBZbCondBandBulkHamiltonian_h_
#define _SBZbCondBandBulkHamiltonian_h_


#include "SBbulkHamiltonian.h"
#include "ZbDDsemiconductor.h"


//! A class to calculate single band Hamiltonian of zinc-blende material
class SBZbCondBandBulkHamiltonian : public SBbulkHamiltonian
{

 public:

  //!Constructor
  SBZbCondBandBulkHamiltonian( );

  //!Constructor
  /*!
    \param parameters parameters that describe a zinc-blende crystal
  */
  SBZbCondBandBulkHamiltonian(ZbDDsemiconductor::ZbDDparameters& parameters);
  
  //! Gamma or X or L
  std::string min_name;

  //! minima number 
  /*!
    O for Gamma; 
    0,1,2 for X;
    0,1,2,3 for L.
  */ 
  unsigned int min_number;

 
  void  calculate_edge_and_mass();

  
  void apply_strain_and_potential(Tensor2Sym& strain_crystal, double el_potential);

 private:

  ZbDDsemiconductor::ZbDDparameters  par; 

};

#endif
