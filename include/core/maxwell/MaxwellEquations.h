#ifndef _MAXWELLEQUATIONS_H_
#define _MAXWELLEQUATIONS_H_

#include "EigenvalueProblem.h"




//!Class to solve Maxwell equations
class MaxwellEquations : public EigenvalueProblem
{
 public:

 
  


  //!constructor
  MaxwellEquations(void) {};
  
  //!destructor
  virtual ~MaxwellEquations(void) {};
  
  virtual PhysicalModel*
    create_physical_model(const ModelOptions& options) const
    throw (ModelErrorException);
    
   
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException) {};

 protected:


 private:

  //!computes Hamiltonian and S matrix
  void calculate_Hamiltonian_and_S(void);

  //!solves eigenvalue problem
  /*!
    \param ev_number number of eigenvalues requested
    \param spectrum_shift additional spetrum shift [eV]
  */
  void solve_eigen_value_problem(unsigned int ev_number, double spectrum_shift = 0.0 );
  

  //! Kronecker delta \f$ \delta_{ij} \f$ 
  inline int delta_Kronecker(int i, int j);
  
  
  //! Tensor product of two Levi-Civita sybmols:  \f$  \sum_{i = 1} ^3 e_{ijk}e_{imn} = \delta_{jm}{kn} - \delta{jn}\delta{km} \f$
  inline int LeviCivita_product(int j, int k , int m, int n);

};



inline int MaxwellEquations::delta_Kronecker(int i, int j)
{
  return ((i==j) ? 1 : 0);
}


//------------------------------------------------------------------------//


inline int MaxwellEquations:: LeviCivita_product(int j, int k , int m, int n)
{
  return ( delta_Kronecker( j,m )*delta_Kronecker( k,n ) -  delta_Kronecker( j,n )*delta_Kronecker( k,m )  );
}

#endif
