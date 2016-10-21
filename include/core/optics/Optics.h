// $Id$

#ifndef _OPTICS_H_
#define _OPTICS_H_

#include "SimulationInterface.h"
#include "EigenvalueProblem.h"
#include "KspaceIntegration.h"

LIBMESHCLASS(Mesh);

class Tensor1;

//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point and compute k-space integrations of the spectrum.
 *
 * k-integrations are performed implementing calculate_for_k_point() 
 * invoked by KspaceIntegrations.
 *
 * In turn this method calls solve_for_kpoint() interface implemented in 
 * the base class EigenvalueProblem for spectral calculations.
 */
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
    options(void) : polariz(0) {}

    double Emin;//<! left boundary of spectrum [eV] 
    double Emax;//<! right boundary of spectrum [eV] 
    double dE;//<! spectrum mesh step [eV]
    double Gamma;//<! spectrum broadening [eV] 
    double nr; //<! the effective refractive index
    Tensor1 polariz; //<! light  polarization [eV] 
    bool  get_occ;  //get occupations from Fermi (Eigenvalues) or set to 1
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
    void calculate_spectrum(const Mesh& Energy, 
		            double Gamma,
			    const Tensor1& polariz, 
                            DofField& spectrum);

    //! Set k-vector for calculation
    void set_k_point(const Point& k_vec);

    //! method invoked by KspaceIntegration via the hook 
    /*!
     * Warning: this method uses internally _spectrum_z as temporary container.
     */
    void calculate_for_k_point(const Point& kpoint, DofField& spectrum, double& estimator);  

    void assemble(const ModelOptions& opt = ModelOptions());

  protected:

    //! init
    virtual void do_init(void);	

    //! solve k.p model
    virtual void do_solve(void);

    //! plot the  spectrum results for the single k-point case
    virtual void do_plot(void); 

    virtual void plot_globaldata();

    virtual void parse_options(void);   

    //! For now do not allow further reimplementation
    void do_setup_solution_variables(void);

    //! For now do not allow further reimplementation
    void get_solution_secure(std::map<ID, std::vector<double> >& values);


    //! assemble the optical matrix problem
    virtual void do_assemble(const ModelOptions& opts) = 0;


    //! compute optical spectrum (currently only for one fixed polarization)
    virtual void do_calculate_spectrum(const Mesh& Energy,
		                       double Gamma,
				       const Tensor1& polariz, 
                                       DofField& spectrum); 


    //! compute matrix elements of P 
    virtual void do_compute_matrix_elements(void) = 0;

    
    //! compute bulk matrix 
    virtual void calculate_matrix_bulk(void) = 0;


    //! numbers of eigensates counting from 0 (e.g. el0, el1, el2, ...)
    std::vector<ID>  _initial_indices;
    
    //! numbers of eigensates counting from 0 (e.g. hl0, hl1, hl2, ...)
    std::vector<ID>  _final_indices; 
    
   
    //! maps the local state indices to EigenstateModel container    
    std::map<ID, ID> _initial_state_numbers;
  
    //! maps the local state indices to EigenstateModel container    
    std::map<ID, ID> _final_state_numbers;     


    //! pointer to the eigenvalue solver for initial states
    EigenvalueProblem* _initial_state_model;
    
    //! pointer to the eigenvalue solver for final states
    EigenvalueProblem* _final_state_model;


    //! type of particle for the intial states 
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
    std::vector <std::vector <libMesh::Complex> >  _P_matrix[3];

    KspaceIntegration* _k_integration;

    //!defines which job has to be done
    JobKind job;

    options _opt; 

  private:

    enum Solutions
    {
      OpticalPower,  //!< the emitted optical power
      Recombination, //!< the total recombination rate
      PeakEnergy     //!< the emission peak energy
    };

    void set_states(void);

    void init_k_space_integration(void);

    void do_k_space_integration(void);

    void compute_matrix_elements(void);

    void print_info(void);

    void check_states(void);

    double lorentzian(double arg, double fwhm);

    double gaussian(double arg, double fwhm);

    //! The total emitted/absorbed power
    double _total_power;

    //! The total recombination/generation rate
    double _recombination;

    //! The emission peak energy
    double _peak_energy;

};



inline
void Optics::calculate_spectrum(const Mesh& Energy, double Gamma,const Tensor1& polariz, 
                            DofField& spectrum)
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


inline
double Optics::lorentzian(double arg, double fwhm)
{
  return (0.5*fwhm / (arg * arg + 0.25 * fwhm * fwhm) / M_PI);
}

inline
double Optics::gaussian(double arg, double fwhm)
{
  double s = fwhm / 2.35482; // 2*sqrt(2*ln(2))
  return (0.5 * M_SQRT1_2 * M_2_SQRTPI / s * exp(-0.5 * arg*arg / (s * s)));
}


#endif // _OPTICS_H_
