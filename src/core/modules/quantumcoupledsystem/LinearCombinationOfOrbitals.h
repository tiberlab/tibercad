/*
 * LinearCombinationOfOrbitals.h
 *
 *  Created on: 18 Jan 2022
 *      Author: miesu
 */

#ifndef _LINEARCOMBINATIONOFORBITALS_H_
#define _LINEARCOMBINATIONOFORBITALS_H_

#include "EigenvalueProblem.h"
#include "EigenSolver.h"



//! Class to calculate electronic properties of large systems of coupled quantum dots-wells

class LinearCombinationOfOrbitals : public EigenvalueProblem
{

  public:

//

    //! Constructor
    LinearCombinationOfOrbitals(const ModelOptions& options);


    //! Destructor
    virtual ~LinearCombinationOfOrbitals(void);


    //! Create object
    static LinearCombinationOfOrbitals* create(const ModelOptions& options);


    //! public member to invoke matrix assembly from other modules
    void assemble(const ModelOptions& options = ModelOptions());



    struct options
      {
        std::string solver;                 //!< solver type

        unsigned int max_iteration_number;  //!< maximum number of iterations for the eigenvalue solver

        double    eigen_solver_tolerance;   //!< tolerance for eigenvalue solver

        unsigned int number_of_eigenstates; //!< number of eigenstates to be calculated

        std::string preconditioner;         //!< preconditioner name

        bool periodicity[3];                //!< periodic boundary conditions

        std::string spectral_trans;         //!< spectral transformation

        double spectrum_shift;              //!< Spectrum shift

        std::string st_ksp_type;            //!< Liner system solution method

        std::string strategy;               //<! matlab (algorithm used in Matlab) or general (recommended by SLEPC)

        bool monitor;   //<! activates convergence monitor if true

        double spectrum_inversion_tolerance; //<! tolerance for spectrum inversion

        bool dump_on_file;

        bool estimate_spectrum_shift; //!< calculate spectrum shift from band edges;

      };



  protected:


    //! init
    virtual void do_init(void) override final;



    virtual int get_H_dim() const override final;



    virtual void get_H_csr(std::vector<libMesh::Complex>& A,
                               std::vector<int>& JA,
                               std::vector<int>& IA,
                               const std::vector<unsigned int>& perm
                               = std::vector<unsigned int>(0)) const override final;



    virtual int get_H_nnz() const override final;


    //! solve k.p model
    virtual void do_solve(void) override final;



    virtual void do_calculate_density_at_k(DofField& Density) override final;



    virtual void do_solve_for_kpoint(const Point& kpoint) override final;


    //assemble final system's Hamiltonian
    virtual void do_assemble(const ModelOptions& options) override final;



    virtual void do_copy_H_to_solver() override final;



    virtual void do_copy_S_to_solver() override final;



    virtual void do_plot(void) override final;



  private:



    EigenvalueProblem* _single_dot_hamiltonian;



    EigenvalueProblem* _total_system_hamiltonian;



    std::vector<size_t> _shift_x;


    //! pointer to the real part of the Hamiltonian
    libMesh::SparseMatrix<double>* _H_lcqo_real;


    //! pointer to the imaginary part of the Hamiltonian
    libMesh::SparseMatrix<double>* _H_lcqo_imag;


    //! pointer to the real part of S matrix
    libMesh::SparseMatrix<double>* _S_lcqo_real;


    //! pointer to the real part of S matrix
    libMesh::SparseMatrix<double>* _S_lcqo_imag;
    
        
    libMesh::NumericVector<double>* total_basis_wave_function_b_real;
    

    libMesh::NumericVector<double>* total_basis_wave_function_b_imag;



    options solver_opt, opt;



    void copy_matrix_to_solver(const char matrix);


    //dimension of the lcqo problem
    unsigned int _H_lcqo_size;


    double k_val; //the k-point for the simulation



    void parse_options(void);



    void solve_eigen_value_problem(unsigned int ev_number);



    bool read_SLEPC_solution(void);



    void postprocess_solution(void);



    void calculate_shift(void);

};



inline
void LinearCombinationOfOrbitals::assemble(const ModelOptions& options)
{
  do_assemble(options);
}




inline
int LinearCombinationOfOrbitals::get_H_dim() const
{
  return 0;
}




inline
int LinearCombinationOfOrbitals::get_H_nnz() const
{
  return 0;
}



inline
void LinearCombinationOfOrbitals::do_copy_S_to_solver()
{
  copy_matrix_to_solver('S');
}



inline
void LinearCombinationOfOrbitals::do_copy_H_to_solver( )
{
  copy_matrix_to_solver('H');
}



inline
void LinearCombinationOfOrbitals::get_H_csr(std::vector<libMesh::Complex>& A,
                               std::vector<int>& JA,
                               std::vector<int>& IA,
                               const std::vector<unsigned int>& perm) const
{
}



inline
void LinearCombinationOfOrbitals::do_plot(void)
{
  _total_system_hamiltonian->plot();
}



#endif // _LINEARCOMBINATIONOFORBITALS_H_
