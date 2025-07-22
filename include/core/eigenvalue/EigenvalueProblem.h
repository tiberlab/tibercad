// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"
#include "KspaceIntegration.h"

#include "libmesh/sparse_matrix.h"
#include <complex>
#include <vector>

namespace libMesh
{
  class UnstructuredMesh;
}

class Kspace;


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
      std::vector<libMesh::Complex> eigen_vector; 
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


    //! Parse e model options for k-space specific stuff
    /*!
     * This method is public, because it is useful for both derived classes
     * and for modules working on top of EigenvalueProblem, like Optics.
     */
    ModelOptions parse_kspace_options(const ModelOptions&);
 
    //! Solve the eigenvalue problem for a certain k-point
    void solve_for_kpoint(const Point& k_point);
    
    //! computes matrix elements between state i of particle_i and state j of particle_j
    virtual libMesh::Complex calculate_matrix_element(const std::string& i_particle,
        unsigned int i, 
        const std::string& j_particle,
        unsigned int j);

    //! Projection on different Brillouin zone
    /*!
     * The main use of this method is for BZ unfolding of supercell dispersions
     * onto the respective primitive cell.
     *
     * \param a the eigenstate to project
     * \param k the k point in the target BZ
     * \return the projection weight
     */
    void project_to_primitive_cell(
        const std::vector<eigen_problem_solution>& a,
        const std::vector<libMesh::Point>& kpoints,
        std::vector<std::vector<libMesh::Complex>>& weights) const;


    //! Project solution states onto specific basis orbitals
    /*!
     * The specific implementation of this method is used to
     * project the solutions onto the different basis functions or orbitals
     * (or atoms, nodes etc) for further insight into the solutions.
     * The projection will be done for each k-point
     *
     * \param bases strings identifying the bases/orbitals to be projected onto
     * \param states the calculated states (solutions)
     * \param projection contains the projections at return
     */
    void project_on_bases(
        const std::vector<std::string>& bases,
        const std::vector<eigen_problem_solution>& states,
        std::vector<std::vector<double>>& projection) const;


    //! get number of states
    unsigned int get_num_states(void) const;

    //! get number of states of a given particle type
    unsigned int get_num_states(const std::string& particle) const;

    //! get states indeces of a given particle type
    std::vector<ID> get_state_indices(const std::string& particle) const;


    //! Write the states to screen
    void write_states(void) const;


    //! passes H matrix to the eigensolver
    void copy_H_to_solver(void);
 
    //! passes S matrix to the eigensolver
    void copy_S_to_solver(void);    
   
    //! get H and S
    virtual int get_H_dim() const = 0;
    
    virtual int get_H_nnz() const = 0;

    //! The Hamiltonian is returned in eV
    /*!
     * The matrix is returned back in CSR format, and is copied.
     * If necessary, a permutation can be provided.
     */
    virtual void get_H_csr(std::vector<libMesh::Complex>& A,
                           std::vector<int>& JA,
                           std::vector<int>& IA,
                           const std::vector<unsigned int>& perm
                                 = std::vector<unsigned int>(0)) const = 0;

    //! The overlap matrix
    /*!
     * Contrary to the Hamiltonian, this must not necessarily
     * be implemented if the problem formulation uses an orthogonal
     * basis. is_generalized() should be used to obtain this
     * information.
     *
     * The matrix is returned back in CSR format, and is copied.
     * If necessary, a permutation can be provided.
     */
    virtual void get_S_csr(std::vector<libMesh::Complex>& A,
                           std::vector<int>& JA,
                           std::vector<int>& IA,
                           const std::vector<unsigned int>& perm
                                 = std::vector<unsigned int>(0)) const;

    virtual void print_H(const std::string& outpath) const;

    //! Return the Hamiltonian Units in eV
    virtual double get_H_units(void) const { return 1.0; }

    //! returns a reference to _solution
    /*!
     * this is dangerous and should be substituted with calls to get_eigenvectors()
     */
    const std::vector<eigen_problem_solution>& get_solution(void) const {return _solution;};

    virtual double get_band_edge(const std::string& band);

    virtual unsigned int get_number_of_bands(void) const { return 0; }

    virtual unsigned int get_degeneracy(void) const { return 1; }


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

    virtual void do_assemble(const ModelOptions& options); 

    virtual void do_plot(void);

    virtual ID do_remember_current_solution(ID id);

    virtual void do_set_to_remembered_solution(ID id);

    virtual void do_delete_remembered_solution(ID id);

    /*!
     * \brief solve the eigenvalue problem
     *
     * \param num_eigenvalues the number of eigenvalues required
     * \param spectrum_shift the spectral shift
     */
    void solve_eigenvalue_problem(unsigned int num_eigenvalues,
                                  double spectrum_shift = 0.0);

    //! read SLEPc solutions
    /*!
     * \return true when all required eigenstates are found
     *
     *  1) Read all eigenvalues
     *  2) Sort the eigenvalues and select those we want
     *  3) Read eigenvectors that correspond to the eigenvalues we want
     *  4) do something else
     *
     *  \return number of eigenvalues to still be computed,
     *      and next spectral shift
    */
    virtual std::pair<unsigned int, double> read_slepc_solution(void);

    //! Implementation of projection on different BZ
    /*!
     * Precondition: \c weights is initialized to the size of
     * \c kpoints, and all elements set to 0.
     */
    virtual void do_project_to_primitive_cell(
        const std::vector<eigen_problem_solution>& a,
        const std::vector<libMesh::Point>& kpoints,
        std::vector<std::vector<libMesh::Complex>>& weights) const;

    //! Implementation of the projection onto specific basis orbitals
    /*!
     * At entry, \c projection container is prepared with the right sizes.
     */
    virtual void do_project_on_bases(
        const std::vector<std::string>& bases,
        const std::vector<eigen_problem_solution>& states,
        std::vector<std::vector<double>>& projection) const;


    //!put spectrum shift energy to be almost equal to the 1st eigenvalue
    virtual double get_new_spectrum_shift(void) { return 0; }

    //! method used to compute quantum dispersion
    virtual void compute_dispersion(const ModelOptions& opts);

    //! method used to plot quantum dispersion
    virtual void plot_dispersion(const std::string& filename);

    //! Write out eigenvalues
    virtual void plot_globaldata(void);

    //! Calculate the DOS
    void calculate_dos(void);

    //! Do k-integration of density
    void integrate_density(DofField& density);

    //! Calculate density
    void calculate_density_at_k(const Point& k_point,
        DofField& density, double& error);

    //! Calculate density for the current solution
    virtual void do_calculate_density_at_k(DofField& density);

    //! Get the k-space
    const Kspace* get_kspace(void) const;


    //! process an element and its neighbours
    void process_element(const Elem* elem, unsigned int entryside,
        std::vector<std::vector<eigen_problem_solution>>& ordered_solutions);

    //! Calculate the scalar product between states a and b
    libMesh::Complex scalar_product(const eigen_problem_solution& a,
                           const eigen_problem_solution& b) const;

    //! Calculate the scalar product between states a and b
    libMesh::Complex scalar_product(const std::vector<libMesh::Complex>& a,
                           const std::vector<libMesh::Complex>& b) const;

    //! pointer to the real part of the Hamiltonian
    libMesh::SparseMatrix<double>* _H_real;

    //! pointer to the imaginary part of the Hamiltonian
    libMesh::SparseMatrix<double>* _H_imag;

    //! pointer to the real part of S matrix 
    libMesh::SparseMatrix<double>* _S_real;

    //! pointer to the real part of S matrix 
    libMesh::SparseMatrix<double>* _S_imag;


    //! Stores the energy values for each k point for dispersions
    /*!
     * The ordering of the k-points is as in k-space mesh (\c _kspace)
     */
    std::vector<std::vector<double>> _dispersion;


    //! Stores projection weights for the k points
    /*!
     * This is filled only in case of unfolding
     */
    std::vector<std::vector<double>> _projection_weights;


    //! The names for the projection onto bases/orbitals/atoms
    std::vector<std::string> _projection_names;

    //! The projections of the different solutions at the different k-points
    /*!
     * This is filled only if projection option is given in input. It will
     * create an additional file, to be used together with the dispersion data.
     * 
     * The first index is the k-point, the second the state index
     */
    std::vector<std::vector<std::vector<double>>> _state_projections;


    //! Indicates whether this is a generalized eigenvalue problem
    bool _haveS;
    
    JobKind _job; //!< a job to do

    //! point for bulk dispersion
    Point _bulk_point;

  private:

    //! A map to contain solutions at k-points
    /*!
     * TODO: this should be transformed into a datastructure that
     * can be automatically serialized to file if it gets large.
     */
    typedef std::map<const Point, std::vector<eigen_problem_solution> > KSolutions;

    bool _new_k;

    libMesh::UnstructuredMesh* _energy_mesh;

    //! Already calculated k-points
    KSolutions _ksolutions;

    //! The reciprocal space description
    Kspace* _kspace;

    //! k-vector in 1/nm
    double _k_vector[3];

    //! to remember solutions
    std::map<ID, std::vector<eigen_problem_solution>> _remembered_sol;


    //! to set that k-vector is not new
    void k_is_old(void);


    //! Reorder states according to a reference solution
    void reorder_states(const std::vector<EigenvalueProblem::eigen_problem_solution>& reference);


    //! Reorder all states to a global order
    /*!
     * States are ordered independent of energy by projecting onto 
     * the nearest k-point. The idea is that coefficients as function
     * of k, c(k), have to be continuous functions of k.
     */
    void reorder_states(const Point& kpoint);

  
};

inline
EigenvalueProblem::EigenvalueProblem(const ModelOptions& options)
 : SimulationInterface(options),
   _energy_mesh(0)
{
  _k_vector[0]=0.0;   _k_vector[1]=0.0;   _k_vector[2]=0.0; 
  _kspace = NULL;
}

inline
void
EigenvalueProblem::assemble(const ModelOptions& options)
{
  do_assemble(options);
}




inline 
void EigenvalueProblem::set_k_point(const Point& k_vec)
{
  for (short i = 0; i < 3; i++) _k_vector[i] = k_vec(i);
  _new_k = true;
}


inline
const Kspace* EigenvalueProblem::get_kspace(void) const
{
  return _kspace;
}


inline
libMesh::Complex
EigenvalueProblem::calculate_matrix_element(const std::string&,
        unsigned int, 
        const std::string&,
        unsigned int)
{
  return 0;
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
bool EigenvalueProblem::is_generalized(void) const
{
  return _haveS;
}

//=======================================================================//

#endif // _EIGENVALUEPROBLEM_H_
