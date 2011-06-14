// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"

#include "sparse_matrix.h"
#include <complex>
#include <vector>
#include <petsc_matrix.h>

//! Abstract class to solve complex valued eigenvalue problem
class EigenvalueProblem : public SimulationInterface
{

  public:

    //! Constructor
    EigenvalueProblem(const ModelOptions& options);

    //! Destructor
    ~EigenvalueProblem(void) { };

    struct eigen_state
    {
      double energy;
      unsigned int index;
    };

    struct eigen_problem_solution
    {
      //! particle type ("electron", "hole", "photon", "exciton", ...) 
      std::string particle;
      //! eigenvalue
      double eigen_energy; 
      //! eigenvector
      std::vector<Complex> eigen_vector; 
      //! statistic type ("Bose" or "Fermi")
      std::string statistics;
      //! electro-chemical potential [eV] \f$ \langle \psi |\mu({\bf r} | \psi \rangle \f$
      double electro_chem_pot;
      //! Level temperature
      double temperature; 
    };

    //! 
    void set_k_vector(std::vector<double> k_point);

    const std::vector<double>& get_k_vector(void) const;
  
    bool has_new_k(void) const;

    void k_is_old(void);


    void get_eigenvalues(const std::string& particle, std::vector<double>& values) const;

    void get_populations(const std::string& particle, std::vector<double>& values) const;

    double get_population(int i) const;

    void assemble(const ModelOptions& options = ModelOptions());
    
    //! computes matrix elements between state i of particle_i and state j of particle_j
    virtual Complex calculate_matrix_element(const std::string& i_particle,
							  unsigned int i, 
							  const std::string& j_particle,
							  unsigned int j){}; 

    //! get number of states of a given particle type
    unsigned int get_num_states(const std::string& particle) const;

    /*! Note: for the moment calculate_matrix_element relays on the fact that the first
     *  n_vb states are for valence, then there are all the electron states.
    */

    void write_states(void) const;

    void write_states(const std::string& filename) const;

    //!passes H matrix to the eigensolver
    void copy_H_to_solver(void);
 
    //!passes S matrix to the eigensolver
    void copy_S_to_solver(void);    
 
    //! returns a reference to _solution
    //! this is dangerous and should be substituted with calls to get_eigenvectors()
    const std::vector<eigen_problem_solution>& get_solution() const
      {return _solution;};

    //! compares eigenstate energy for electrons 
    static bool compare_eigen_energy_electrons(const eigen_state& state1, const eigen_state& state2);

    //! compares eigenstate energy for holes 
    static bool compare_eigen_energy_holes(const eigen_state& state1, const eigen_state& state2);
   
  protected:

    std::vector<eigen_problem_solution> _solution;

    double Fermi(double Energy, double Fermi_energy, double Temperature) const;

    double Bose(double Energy, double elec_chem, double Temperature) const;


    virtual void do_copy_H_to_solver(void){};

    virtual void do_copy_S_to_solver(void){};  

    virtual void do_assemble(const ModelOptions& options){};  

    //!read SLEPc solutions
    /*!
 
    1) Read all eigenvalues
    2) Sort the eigenvalues and select those we want 
    3) Read eigenvectors that correspond to the eigenvalues we want
    4) do something else
 
    \param number_of_ev number of eigen functions to read
    */
    virtual void read_SLEPC_solution(unsigned int number_of_ev) {};


    //!put spectrum shift energy to be almost equal to the 1st eigenvalue
    virtual double get_new_spectrum_shift(void){};


    //!pointer to the real part of the Hamiltonian
    //SparseMatrix<Complex>* _H;

    //!pointer to the imaginary part of the Hamiltonian
    SparseMatrix<double>* _H_real;

    //!pointer to the imaginary part of the Hamiltonian
    SparseMatrix<double>* _H_imag;

    //!pointer to the real part of S matrix 
    //SparseMatrix<Complex>* _S;

    //!pointer to the real part of S matrix 
    SparseMatrix<double>* _S_real;

    //!pointer to the real part of S matrix 
    SparseMatrix<double>* _S_imag;

  private:

  std::vector<double> _k_vector;
  
  bool _new_k;

};


inline
EigenvalueProblem::EigenvalueProblem(const ModelOptions& options)
 : SimulationInterface(options)
{
  _k_vector.reserve(3);
  _k_vector[0]=0.0;   _k_vector[1]=0.0;   _k_vector[2]=0.0; 
}

inline
void
EigenvalueProblem::assemble(const ModelOptions& options)
{
  do_assemble(options);
}

inline
const std::vector<double>& EigenvalueProblem::get_k_vector() const
{
  return _k_vector;
}

inline 
void EigenvalueProblem::set_k_vector(std::vector<double> k_point)
{
  _k_vector = k_point;
  _new_k = true;
}

inline 
bool EigenvalueProblem::has_new_k(void) const
{
  return _new_k;
}

inline
void EigenvalueProblem::k_is_old(void)
{
  _new_k=false;
}

inline
bool EigenvalueProblem::compare_eigen_energy_holes(const eigen_state& state1, const eigen_state& state2)
{
  return(state1.energy> state2.energy);
}

//=======================================================================//

inline
bool EigenvalueProblem::compare_eigen_energy_electrons(const eigen_state& state1, const eigen_state& state2)
{
  return(state1.energy< state2.energy);
}


//=======================================================================//

#endif // _EIGENVALUEPROBLEM_H_
