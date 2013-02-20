// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"
#include "Kspace.h"

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

  
    //! Eigenstate structure useful to the solver when sorting states
    struct eigen_state
    {
      double energy;
      unsigned int index;
    };

    //! Container for eigenstates 
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

    //! Set k-vector for calculation
    void set_k_point(const Point& k_vec);

    //! to check if k-vector has changed 
    bool has_new_k(void) const;

    //! to set that k-vector is not new
    void k_is_old(void);

    //! Used to retrieve eigenvalues from other modules
    void get_eigenvalues(const std::string& particle, std::vector<double>& values) const;

    //! Used to retrieve populations from other modules
    void get_populations(const std::string& particle, std::vector<double>& values) const;

    //! get population of a single state   
    double get_population(int i) const;

    //! public member to invoke matrix assembly from other modules
    void assemble(const ModelOptions& options = ModelOptions());
 
    void solve_for_kpoint(const Point& k_point);
    
    //! computes matrix elements between state i of particle_i and state j of particle_j
    virtual Complex calculate_matrix_element(const std::string& i_particle,
							  unsigned int i, 
							  const std::string& j_particle,
							  unsigned int j){}; 

    //! get number of states
    unsigned int get_num_states(void) const;

    //! get number of states of a given particle type
    unsigned int get_num_states(const std::string& particle) const;

    //! get states indeces of a given particle type
    std::vector<unsigned int> get_state_indices(const std::string& particle) const;


    /* Note: for the moment calculate_matrix_element relays on the fact that the first
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
    const std::vector<eigen_problem_solution>& get_solution(void) const {return _solution;};

    //! compares eigenstate energy for electrons needed for sorting
    static bool compare_eigen_energy_electrons(const eigen_state& state1, const eigen_state& state2);

    //! compares eigenstate energy for holes needed for sorting
    static bool compare_eigen_energy_holes(const eigen_state& state1, const eigen_state& state2);
   
  protected:

    std::vector<eigen_problem_solution> _solution;

    double Fermi(double Energy, double Fermi_energy, double Temperature) const;

    double Bose(double Energy, double elec_chem, double Temperature) const;

    virtual void init_kspace(void);

    virtual void do_solve_for_kpoint(const Point& k_point){};

    virtual void do_copy_H_to_solver(void){};

    virtual void do_copy_S_to_solver(void){};  

    virtual void do_assemble(const ModelOptions& options){}; 

    virtual void do_plot(void);

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

    ModelOptions parse_kspace_options(const ModelOptions&);

    //! method used to compute quantum dispersion
    virtual void compute_dispersion(void);

    //! method used to plot quantum dispersion
    virtual void plot_dispersion(void);

    //!pointer to the imaginary part of the Hamiltonian
    SparseMatrix<double>* _H_real;

    //!pointer to the imaginary part of the Hamiltonian
    SparseMatrix<double>* _H_imag;

    //!pointer to the real part of S matrix 
    SparseMatrix<double>* _S_real;

    //!pointer to the real part of S matrix 
    SparseMatrix<double>* _S_imag;

    //!k-vector in atomic units
    double _k_vector[3];

    Kspace* _kspace;

    bool do_dispersion;

    std::vector< std::vector<double> > _dispersion;

    int disp_range[2];

  private:
  
     bool _new_k;

};

inline
EigenvalueProblem::EigenvalueProblem(const ModelOptions& options)
 : SimulationInterface(options)
{
  _k_vector[0]=0.0;   _k_vector[1]=0.0;   _k_vector[2]=0.0; 
  _kspace = NULL;
  do_dispersion = false;
  disp_range[0]=0; 
  disp_range[1]=0;
}

inline
void
EigenvalueProblem::assemble(const ModelOptions& options)
{
  do_assemble(options);
}

//inline
//std::vector<eigen_problem_solution>& EigenvalueProblem::get_solution(void) const 
//{
//   return _solution;
//}

inline 
void EigenvalueProblem::set_k_point(const Point& k_vec)
{
  for (short i = 0; i < 3; i++) _k_vector[i] = k_vec(i);
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
