// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"
#include "KspaceIntegration.h"

#include "sparse_matrix.h"
#include <complex>
#include <vector>


class Kspace;
class Mesh;


//! Abstract class to solve complex valued eigenvalue problem
class EigenvalueProblem : public SimulationInterface
{

  public:

    //!control options
    enum JobKind
    {
      EIGENSTATES = 0, //!< eigenenergies
      BULKEIGENSTATES=1, //!< bulk eigenenergies
    };


    //! Constructor
    EigenvalueProblem(const ModelOptions& options);

    //! Destructor
    ~EigenvalueProblem(void) { };

  
    //! Eigenstate structure useful to the solver when sorting states
    struct eigen_state
    {
      double energy;
      unsigned int index;
      std::string particle;
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

    //! Get the k point
    Point get_k_point(bool relative_coord = false) const;

    //! to check if k-vector has changed 
    bool has_new_k(void) const;

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
							  unsigned int j) { return 0; }

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


    //!passes H matrix to the eigensolver
    void copy_H_to_solver(void);
 
    //!passes S matrix to the eigensolver
    void copy_S_to_solver(void);    
   
    //! get H and S
    virtual int get_H_dim() const { return 0; }
    
    virtual int get_H_nnz() const { return 0; }

    //! The Hamiltonian is returned in eV 
    virtual void get_H_csr(std::vector<Complex>& A,std::vector<int>& JA,std::vector<int>& IA) const {};

    virtual void get_S_csr(std::vector<Complex>& A, std::vector<int>& JA,std::vector<int>& IA) const {};

    virtual void print_H(const std::string& outpath) const {};

    //! Return the Hamiltonian Units in eV
    virtual double get_H_units(void) const{};

    //! returns a reference to _solution
    //! this is dangerous and should be substituted with calls to get_eigenvectors()
    const std::vector<eigen_problem_solution>& get_solution(void) const {return _solution;};

    //  ! compares eigenstate energy for electrons needed for sorting
    //static bool compare_eigen_energy_electrons(const eigen_state& state1, const eigen_state& state2);

    //  ! compares eigenstate energy for holes needed for sorting
    //static bool compare_eigen_energy_holes(const eigen_state& state1, const eigen_state& state2);

    virtual double get_band_edge(const std::string& band) { return 0; }

    virtual unsigned int get_number_of_bands(void) const { return 0; }

    void set_permutation(const std::vector<unsigned int>& p);

    void init_permutation(const unsigned int n_dofs);

    bool is_generalized() const;

    /*!
     * \brief initialize the solution container to hold \c num_solutions solutions
     *
     * It is responsibility of the user to call this method before solving
     * or before extracting solutions.
     */
    void initialize_solution_container(size_t num_solutions);

     
  protected:

    std::vector<eigen_problem_solution> _solution;

    double Fermi(double Energy, double Fermi_energy, double Temperature) const;

    double Bose(double Energy, double elec_chem, double Temperature) const;

    virtual void init_kspace(const ModelOptions& opt);

    virtual void do_solve_for_kpoint(const Point& k_point);

    virtual void do_copy_H_to_solver(void){};

    virtual void do_copy_S_to_solver(void){};  

    virtual void do_assemble(const ModelOptions& options){}; 

    virtual void do_plot(void);

    virtual ID do_remember_current_solution(ID id);

    virtual void do_set_to_remembered_solution(ID id);

    virtual void do_delete_remembered_solution(ID id);

    //!read SLEPc solutions
    /*!
     * \return true when all required eigenstates are found
 
    1) Read all eigenvalues
    2) Sort the eigenvalues and select those we want 
    3) Read eigenvectors that correspond to the eigenvalues we want
    4) do something else
 
    \param number_of_ev number of eigen functions to read
    */
    virtual bool read_SLEPC_solution(void) { return true; }


    //!put spectrum shift energy to be almost equal to the 1st eigenvalue
    virtual double get_new_spectrum_shift(void) { return 0; }

    ModelOptions parse_kspace_options(const ModelOptions&);

    //! method used to compute quantum dispersion
    virtual void compute_dispersion(void);

    //! method used to plot quantum dispersion
    virtual void plot_dispersion(const std::string& filename);

    //! Write out eigenvalues
    virtual void plot_globaldata(void);

    //! Calculate the DOS
    void calculate_dos(void);

    //! process an element and its neighbours
    void process_element(const Elem* elem, unsigned int entryside,
        std::vector<std::vector<eigen_problem_solution>>& ordered_solutions);

    //! Calculate the scalar product between states a and b
    Complex scalar_product(const eigen_problem_solution& a,
                           const eigen_problem_solution& b) const;

    //! Calculate the scalar product between states a and b
    Complex scalar_product(const std::vector<Complex>& a,
                           const std::vector<Complex>& b) const;

    //! pointer to the real part of the Hamiltonian
    SparseMatrix<double>* _H_real;

    //! pointer to the imaginary part of the Hamiltonian
    SparseMatrix<double>* _H_imag;

    //! pointer to the real part of S matrix 
    SparseMatrix<double>* _S_real;

    //! pointer to the real part of S matrix 
    SparseMatrix<double>* _S_imag;

    //! Stores a general permutation on dofs
    std::vector<ID> _perm;
    std::vector<ID> _inv_perm;

    //! bool do_dispersion;
    std::vector< std::vector<double> > _dispersion;

    //int disp_range[2];
    bool _haveS;
    
    JobKind _job; //!< a job to do

    //!point for bulk dispersion
    Point _bulk_point;

  private:

    //! A map to contain solutions at k-points
    typedef std::map<const Point, std::vector<eigen_problem_solution> > KSolutions;

    bool _new_k;

    Mesh* _energy_mesh;

    //! Already calculated k-points
    KSolutions _ksolutions;

    //! The reciprocal space description
    Kspace* _kspace;

    //! k-vector in arbitrary units (for now)
    double _k_vector[3];

    //! to remember solutions
    std::map<ID, std::vector<eigen_problem_solution>> _remembered_sol;

    //! Callback function for dos_calculation
    void _dos_for_kpoint(const Point& k_point,
        const Point& refpoint,
        DofField& density,
        double& integrated_quantity);


    //! to set that k-vector is not new
    void k_is_old(void);

  
};

inline
EigenvalueProblem::EigenvalueProblem(const ModelOptions& options)
 : SimulationInterface(options),
   _energy_mesh(0)
{
  _k_vector[0]=0.0;   _k_vector[1]=0.0;   _k_vector[2]=0.0; 
  _kspace = NULL;
  //do_dispersion = false;
  //disp_range[0]=0;
  //disp_range[1]=0;
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
void EigenvalueProblem::init_permutation(const unsigned int n_dofs)
{
  _perm.resize(n_dofs);
  for(unsigned int i=0; i<n_dofs; i++)
    _perm[i] = i;
}

inline
void EigenvalueProblem::set_permutation(const std::vector<unsigned int>& p)
{
  _perm.resize(p.size());
  for(unsigned int i=0; i<p.size(); i++)
    _perm[i] = p[i];
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

//inline
//bool EigenvalueProblem::compare_eigen_energy_holes(const eigen_state& state1, const eigen_state& state2)
//{
//  return(state1.energy> state2.energy);
//}

//=======================================================================//

//inline
//bool EigenvalueProblem::compare_eigen_energy_electrons(const eigen_state& state1, const eigen_state& state2)
//{
//  return(state1.energy< state2.energy);
//}

inline
bool EigenvalueProblem::is_generalized(void) const
{
  return _haveS;
}

//=======================================================================//

#endif // _EIGENVALUEPROBLEM_H_
