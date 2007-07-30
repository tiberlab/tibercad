#ifndef _MAXWELLEQUATIONS_H_
#define _MAXWELLEQUATIONS_H_



class MaxwellEquations : public SimulationInterface
{
 public:

 
  


  //!constructor
  MaxwellEquations(void);
  
  //!destructor
  virtual ~MaxwellEquations(void);
  
  

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
  


};



#endif
