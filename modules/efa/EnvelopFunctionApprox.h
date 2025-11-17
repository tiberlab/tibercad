// $Id$

#ifndef _ENVELOPFUNCTIONAPPROX_H_ 
#define _ENVELOPFUNCTIONAPPROX_H_


#include "FEMEigenvalueProblem.h"
#include "TemperatureInterface.h"
#include "StrainInterface.h"


//! A class that constructs EFA Hamiltonian and S-matrix 
class EnvelopFunctionApprox  : public FEMEigenvalueProblem
{
 public:

  //! data structure that contains options for effective mass
  struct options
  {
    
    unsigned int num_el_states; //!< number of electron states to be calculated

    unsigned int num_hl_states; //!< number of hole states to be calculated

    unsigned int degeneracy; //!< the degeneracy factor

    bool  consider_potential; //!< apply potential to the EFA Hamiltonian;

    bool  consider_potential_bulk; //!< apply potential to the bulk EFA Hamiltonian;

    bool  consider_strain_bulk; //!< apply strain to the bulk EFA Hamiltonian;

    bool estimate_spectrum_shift; //!< calculate spectrum shift from band edges;

    bool convergent_density;//!< if true, the number of eigenstates will be increased to reach the tolerance
 
    unsigned int initial_eigenstates_number; //!< initial number of eigenstates that is used in an iterative calculation of the density

    double relative_density_tolerance; //!< stops itarations if \f$ \rho_i / \rho_{i+1} < \varepsilon    \f$, where \f$ \rho \f$ is the                                              total density 
 
    double eigen_number_increase_factor; //!< to increase number of eigenstates for the next iteration 

    bool local_occupation; //!<If a local occupation is considered 

    unsigned int first_state;

    double k_val;
    
    bool assume_paraboloid;

    bool analytic_k_int; //!< do analytic k-integration, assuming parabolic DOS


  };



 
  //!constructor
  EnvelopFunctionApprox(const ModelOptions& options);


  //!destructor
  ~EnvelopFunctionApprox();

 
 
  //! returns conduction band minima for electrons and valence band maximum for holes 
  double get_band_edge(const std::string& band);


  
  //!calculate averaged value of the electrochemical potential \f$\langle \psi|\mu|psi \rangle \$f
  /*!
    \param  i number of state
   
   */
  double calculate_fermi_averaged(unsigned int i, const std::string& particle);


  
  //!calculate averaged value of the temperature \f$\langle \psi|T|psi \rangle \$f
  /*!
    \param  i number of state
   
   */
  double calculate_temperature_averaged(unsigned int i);
  
 


  //! sets opt.initial_eigestates_number
  void set_initial_eigenstates_number(unsigned int n);



  //! returns number of active cells
  unsigned int get_number_of_active_cells();


  static  EnvelopFunctionApprox* create(const ModelOptions& options);


  virtual PhysicalModel*
    create_bulk_model(const ModelOptions& options,
        const Material* mat) const;
    
   

  //!returns number of bands
  unsigned int get_number_of_bands(void) const;

  unsigned int get_degeneracy(void) const;

  //!returns map for kp bands
  inline const std::map<short, short>& get_band_map(void) const;


  //! get the  vector of  eigenvalues 
  void get_eigenenergies(std::vector<double>& values) const;

  //! get occupation of eigenstates
  void get_occupations(std::vector<double>& values) const;


  //! returns element used for bulk calculation 
  inline const Elem*  return_bulk_element(void) const;

 private:

  enum Variables
  {
    ProbabilityDensity,
    EigenEnergy,
    Occupation,
    EigenEnergyOnMesh,
    eDensity,
    hDensity,
    EnvelopeFunctions
  };

  //! interface for temperature acquisition
  TemperatureInterface _temp_interface;
 
  //! Interface to strain simulation
  SolutionProvider _strain_interface;

  options opt;

  //!< number of bands in EFA 
  unsigned int number_of_bands; 

  //!< map between band numbers: from 8 band scheme to any band scheme
  std::map<short, short> band_map; 
  
  std::string system_name;

  //! The solution provider for the electrostatic potential
  SolutionProvider _el_pot;

  //! The solution provider for the electron electrochemical potential
  SolutionProvider _el_elchem;

  //! The solution provider for the hole electrochemical potential
  SolutionProvider _hl_elchem;

  //! The solution provider for the CB band edge
  SolutionProvider _cb_edge;

  //! The solution provider for the VB band edge
  SolutionProvider _vb_edge;

  //!pointer to a drift-diffusion object that is used to get potential data 
  SimulationInterface* poisson_equation;


  //! ID of electric potential in poisson_equation
  ID potential_ID;

  //! ID of electro chemical potential in poisson_equation for particles
  ID el_electro_chem_pot_ID;
  ID hl_electro_chem_pot_ID;

  //! band edge ID in poisson_equation
  ID cb_band_edge_ID;
  ID vb_band_edge_ID;

  
  //! If \c true we have to calculate the density
  bool _calculate_density;
  

  //! The quadrature rule
  libMeshEnums::QuadratureType _quadrature_type;

  //! The square root of the inverse of the overlap matrix
  /*!
   * This is only used for trapezoidal rule that produces a diagonal
   * overlap matrix.
   */
  std::vector<double> _sqrt_S_inv;


  //!bands names
  std::vector<std::string> psi_name;
 

  //!calculates analytical charge density for a 1D/2D structures
  /*!
    \f$  \rho({\bf r}) = |\psi({\bf r})|^2 \frac{mkT}{2 \pi \hbar^2}\ln (1 + \exp (\frac{\mu - E}{kT}) )    \f$
    \f$  \rho({\bf r}) = |\psi({\bf r})|^2 \sqrt{\left( \frac{mkT}{2\pi\hbar^2} \right)}
    F_{-1/2}\left(  \frac{\mu - E}{kT}         \right)       \f$
  */
  void calculate_density_analytic(void);


  //!read SLEPc solutions
  /*!
 
  1) Read all eigenvalues
  2) Sort the eigenvalues and select those we want 
  3) Read eigenvectors that correspond to the eigenvalues we want
  4) normalize eigenfunctions
  5) calculate fermi energy for each state
 
    return number of eigenvalues to still be computed,
           and next spectral shift
  */
  std::pair<unsigned int, double> read_slepc_solution(void) override;


  //! Transform the eigenstate with S^-1/2, if needed
  void transform_eigenstate(std::vector<libMesh::Complex>& eigvec);

 
  //!number of nodes used in the model
  unsigned int number_of_nodes;


  //! compares eigenstate energy for electrons 
  static bool compare_eigen_energy_electrons(const double& state1, const double& state2);

  //! compares eigenstate energy for holes 
  static bool compare_eigen_energy_holes(const double& state1, const double& state2);
  

  //!list of periodic nodes
  std :: vector< std :: vector <const Node*> >  nodes_periodic; //dim node list's: each contains list of nodes that periodic b.c
                                                                //must be applied to
  

  void calculate_Hamiltonian_and_S(void);


  //! calculates the norm of the eigenstate \f$ \sqrt {| \langle \psi|\psi \rangle |} \f$
  /*!
    \param state_number number of the eigenstate
    \param projections the projectiosn on the single bloch states, if required
  */
  double eigenstate_norm(unsigned int state_number, std::vector<double> *projections = nullptr) const;


  //!Calculates Fermi Dirac probability
  /*!
    For electrons:  \f$ p = \frac{1}{1 + \exp (\frac{E - \mu}{kT})}     \f$
    For holes:      \f$ p = 1 - \frac{1}{1 + \exp (\frac{E - \mu}{kT})} \f$
    \param Energy   state energy [eV]
    \param Fermi_energy  Fermi energy [eV]
    \param temperature   temperature [K] 
   */
  double Fermi_statistics_probability(double Energy, double Fermi_energy,
      double temperature, const std::string& particle) const;


  //!Calculate number of bands in the Hamiltonian 
  short calculate_number_of_bands(void) const;



  //!calculated density
  std::map<const Elem*, double> _density;
 
   //! element used for bulk calculations
  const Elem* _bulk_mat_element;

  //!estimates shift according to band-edge
  void estimate_spectrum_shift(void);
  

  //!put spectrum shift energy to be almost equal to the 1st eigenvalue
  double get_new_spectrum_shift(void);


  //!return electric potential (wrapper)
  void  get_electric_potential(const Elem* elem, const std::vector<Point>& q_point, std::vector<double> electric_potential) const;

  //!returns band edge (without potential)
  double get_band_edge(const Elem* elem, const std::string& particle) const;

  //!return electric potential from drift-diffusion
  double get_electric_potential(const Elem* elem, const Point&  qp) const;

  //!return electro chemical potential from drift-diffusion
  double get_el_electro_chem_potential(const Elem* elem) const;
  double get_hl_electro_chem_potential(const Elem* elem) const;




  //! bulk eigenstates
  void solve_bulk(void);

  virtual double get_H_units(void) const final;

  //! Redeclare the solutions, if number of states changed
  void redeclare_solutions(void);

 protected:


  virtual void get_solution_secure(const Elem* elem,
      std::map<ID, std::vector<double> >& values,
      const std::vector<Point>& points) final;

  virtual void get_solution_secure(
      std::map<ID, std::vector<double> >& values) final;


  virtual void do_setup_solution_variables(void) final;

  //! set to a remembered solution
  /*
   * The actual implementation is in a base class, but we
   * need to also redeclare the solution variables with the
   * right number of solutions.
   */
  virtual void do_set_to_remembered_solution(ID id) final;

  virtual void 	do_init(void) final;

  virtual void 	do_solve (void) final;

  virtual void do_print_info(void) final;

  virtual void 	parse_options (void) final;

  //virtual void do_plot(void);

  virtual void do_calculate_density_at_k(DofField& density) final;


  virtual void do_solve_for_kpoint(const Point& kpoint) final;

 
  virtual void do_assemble(const ModelOptions& options) final;


  //! Implements projection onto Bloch basis functions
  virtual void do_project_on_bases(
      const std::vector<std::string>& bases,
      const std::vector<eigen_problem_solution>& states,
      std::vector<std::vector<double>>& projection) const final;



 private:

  bool check_confinement(const std::vector<libMesh::Complex>& state);

  double eigenstate_norm(const std::vector<libMesh::Complex>& eigen_vector,
                         std::vector<double> *projections = nullptr) const;

  //! Get strain in crystal coordinates
  void get_crystal_strain(const libMesh::Elem* elem,
      const libMesh::Point& point, Tensor2& strain);

};


//-------------------------------------------------------------------
inline double EnvelopFunctionApprox::Fermi_statistics_probability(double Energy, double Fermi_energy,
                                                                  double Temperature, const std::string& particle) const
{
  

  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  if (particle != "el")
    exp_arg *= -1.0;

  double el_fermi;

  if (exp_arg > 20)
    el_fermi = std::exp(-exp_arg);
  else
    el_fermi = 1.0/(  1.0 +  std::exp(exp_arg)  );

  return(el_fermi);
}

//-------------------------------------------------------------------

//---------------------------------------------------------

inline void EnvelopFunctionApprox::set_initial_eigenstates_number(unsigned int n)
{
  opt.initial_eigenstates_number = n;
}


//---------------------------------------------------------

inline EnvelopFunctionApprox*  EnvelopFunctionApprox::create(const ModelOptions& options)
{
  return (new EnvelopFunctionApprox(options) );
}


//---------------------------------------------------------

inline unsigned int EnvelopFunctionApprox::get_number_of_bands(void) const
{
  return (number_of_bands);
}

//---------------------------------------------------------

inline unsigned int EnvelopFunctionApprox::get_degeneracy(void) const
{
  return (opt.degeneracy);
}

//--------------------------------------------------------------------
inline const std::map<short, short>& EnvelopFunctionApprox::get_band_map(void) const
{

  return (band_map);

}



inline const Elem*  EnvelopFunctionApprox::return_bulk_element(void) const
{
  return _bulk_mat_element;
}

inline double EnvelopFunctionApprox::get_H_units(void) const
{
  return Constants::Hartree;
}

#endif
