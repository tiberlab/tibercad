// $Id$

#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "SimulationInterface.h"
#include "EigenvalueProblem.h"

class Mesh;
class Elem;
class Tensor1;
class KspaceIntegration;
//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point and compute k-space integrations of the spectrum
 */
typedef std::map<const Elem*, double> DofField;


class Optics : public SimulationInterface
{

  public:

  //!control options
  enum JobKind
  {
    MATREL = 0, //!< calculate matrix elements of momentum matrix
    BULKMATREL = 1, //!<  calculate matrix elements of momentum matrix for bulk simulation
  };

  struct options 
  {
    double Emin;//<! left boundary of spectrum [eV] 
    double Emax;//<! right boundary of spectrum [eV] 
    double dE;//<! spectrum mesh step [eV]
    double Gamma;//<! spectrum broadening [eV] 
    Tensor1 polariz; //<! light  polarization [eV] 
  };


    //! The constructor
    Optics(const ModelOptions& options);

    //! The destructor
    virtual ~Optics(void);


    //!calculate spectrum 
    /*!
     \f$
    
     P(\hbar \omega) = \sum_{i,j} \frac{1}{2\pi^2}  \frac{\omega^2_{ij} e^2 }{m^2 c^3}  
             |{\bf M_{i,j} e}|^2 f_i(E_i)(1 - f_j(E_j)) 
     \frac{\Gamma/2} {(\hbar \omega_{ij} - \hbar \omega)^2 + (\Gamma/2)^2} d\Omega
     \f$
      
     \param Energy energy grid [eV]
     \param spectrum calculated spectrum (atomic units)
     \param Gamma broadering parameter [eV]
     \param polariz polarization vector of a linearly polarized light (must be a normalized one, \f$ |{\bf e}| = 1 \f$)
    
    */
    void calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
                                                std::map<const Elem*, double>& spectrum);

    //! Set k-vector for calculation
    void set_k_point(const Point& k_vec);

    //! call-back method that KspaceIntegration invokes
    void calculate_for_k_point(const Point& kpoint, DofField& spectrum, double& estimator);  

    void assemble(const ModelOptions& opt = ModelOptions());

  protected:

    //! init
    virtual void do_init(void);	

    //! solve k.p model
    virtual void do_solve(void);

    //! plot the  spectrum results for the single k-point case
    virtual void do_plot(void); 


    virtual void parse_options(void);   


    //! assemble the optical matrix problem
    virtual void do_assemble(const ModelOptions& opts) = 0;


    //! compute optical spectrum (currently only for one fixed polarization)
    virtual void do_calculate_spectrum(const Mesh& Energy, double Gamma, const Tensor1& polariz, 
                                    std::map<const Elem*, double>& spectrum); 


    //! compute matrix elements of P 
    virtual void do_compute_matrix_elements(void) = 0;

    
    //! compute bulk matrix 
    virtual void calculate_matrix_bulk(void) = 0;


    //!numbers of eigensates that are considered as intial states for optical transition
    std::vector<unsigned int> _initial_eigen_state_numbers;
    
    
    //!numbers of eigensates that are considered as final states for optical transition
    std::vector<unsigned int> _final_eigen_state_numbers;
    

    //!pointer to the eigenvalue solver for initial states
    EigenvalueProblem* _initial_state_model;
    
    
    //!pointer to the eigenvalue solver for final states
    EigenvalueProblem* _final_state_model;


    //!type of particle for the intial states 
    std::string _initial_state_particle;


    //!type of particle for the final states 
    std::string _final_state_particle;
    
    //! initial states
    std::vector<EigenvalueProblem::eigen_problem_solution> _i_states;    

    //! final states
    std::vector<EigenvalueProblem::eigen_problem_solution> _f_states;  

    //!k-vector in atomic units
    double _k_vector[3];
 
    Mesh* _energy_mesh;

    DofField _spectrum_x;

    DofField _spectrum_y;

    DofField _spectrum_z;

    //!  momentum matrix elements (should go in Optics.h)
    /*!
      _P_matrix[p][i][j]:  i - initial state; j - final state
    */    
    std::vector <std::vector <Complex> >  _P_matrix[3];

    KspaceIntegration* _k_integration;

    //!defines which job has to be done
    JobKind job;

    options _opt; 

  private:

    void set_states(void);

    void init_k_space_integration(void);

    void do_k_space_integration(void);

    void compute_matrix_elements(void);

    void print_info(void);

    void check_states(void);

};



inline
void Optics::calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
                            std::map<const Elem*, double>& spectrum)
{
  do_calculate_spectrum(Energy, Gamma, polariz, spectrum);
}

inline 
void Optics::set_k_point(const Point& k_vec)
{
  for (short i = 0; i < 3; i++) _k_vector[i] = k_vec(i);
}

inline
void Optics::assemble(const ModelOptions& opts)
{ 
   do_assemble(opts);
}

#endif // _OPTICS_H_
