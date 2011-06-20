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
  //!control options

  enum JobKind
  {
    EIGENSTATES = 0, //!< eigenenergies
    DENSITY = 1, //!<particle density
    BULKEIGENSTATES=2, //!< bulk eigenenergies
    BULKDENSITY = 3 //!< bulk density
  };




  //! data structure that contains options for effective mass
  struct options
  {
    
    unsigned int number_of_bands; //!< number of bands in EFA 

    std::map<short, short> kp_bands; //!< map between band numbers: from 8 band scheem to any band scheem

    std::string particle;   //!< particle name "el" or "hl"

    unsigned int degeneracy; //!< the degeneracy factor

    //double  length_scale;   //!< mesh length scale [Bohr radius]

    //bool periodicity[3];    //!< periodic boundary conditions

    double spectrum_shift;    //!< shift of spectrum ised in matrix assembly[eV]

    bool  consider_potential; //!< apply strain effect to the EFA Hamiltonian;

    bool estimate_spectrum_shift; //!< calculate spectrum shift from band edges;

    bool convergent_density;//!< if true, the number of eigenstates will be increased to reach the tolerance
 
    unsigned int initial_eigenstates_number; //!< initial number of eigenstates that is used in an iterative calculation of the density

    double relative_density_tolerance; //!< stops itarations if \f$ \rho_i / \rho_{i+1} < \varepsilon    \f$, where \f$ \rho \f$ is the                                              total density 
 
    double eigen_number_increase_factor; //!< to increase number of eigenstates for the next iteration 

    JobKind job; //!< a job to do

    bool local_occupation; //!<If a local occupation is considered 

    unsigned int first_state;

    double k_val;
    
    bool assume_paraboloid;


  };



  //struct eigen_energy
  //{
  //  double energy; //!< eigen energy [eV]
  //  unsigned int global_number; //< eigen vector
  //};



  //! data structure that contains eigenvalue and eigenvector
  //struct eigen_propblem_solution
  //{
  //  double eigen_energy; //!< eigen energy [eV]
  //  std::vector< std::complex<double>  > eigen_vector; //< eigen vector
  //  double Fermi_energy; //< electro-chemical potential [eV] \f$ \langle \psi |\mu({\bf r} | \psi \rangle  \f$
  //  double Temperature; //!< averaged temperature for the state [K] 
  //};

  //! returns a reference to solutions
  //const std::vector<eigen_propblem_solution>& get_solution() const;
 
  //!constructor
  EnvelopFunctionApprox(const ModelOptions& options);


  //!destructor
  ~EnvelopFunctionApprox();

 
 
  //! returns conduction band minima for holes and valence band maximum for holes 
  double get_band_edge() const;


  
  //!calculate averaged value of the electrochemical potential \f$\langle \psi|\mu|psi \rangle \$f
  /*!
    \param  i number of state
   
   */
  double calculate_fermi_averaged(unsigned int i);


  
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
    
   

  //!returns a constant reference calculated density
  const std::map<const Elem*, double>& get_density(void) const 
  { return _density; } ;

  //!returns number of bands
  inline unsigned int get_number_of_bands(void) const;


  //!returns map for kp bands
  inline const std::map<short, short>& get_kp_bands(void) const;


  //! get the  vector of  eigenvalues 
  void get_eigenenergies(std::vector<double>& values) const;

  //! get occupation of eigenstates
  void get_occupations(std::vector<double>& values) const;


  inline double get_particle_charge(void) const;
 
  //! element used for bulk calculations
  const Elem* _bulk_mat_element;

  //! returns element used for bulk calculation 
  inline const Elem*  return_bulk_element(void) const;

 private:

  enum Variables
  {
    ProbabilityDensity,
    EigenEnergy,
    Occupation,
    EigenEnergyOnMesh,
    QuantumDensity
  };

  //! interface for temperature acquisition
  TemperatureInterface _temp_interface;
 
  //! Interface to strain simulation
  StrainInterface _strain_interface;

  options opt;

  
  std::string system_name;

  //!pointer to a drift-diffusion object that is used to get potential data 
  SimulationInterface* poisson_equation;


  //! ID of electric potential in poisson_equation
  ID potential_ID;

  //! ID of electro chemical potential in poisson_equation for particles
  ID electro_chem_pot_ID;

  //! band edge ID in poisson_equation
  ID band_edge_ID;

  
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
 

  void set_k_vector(const RealVectorValue& k_vec);

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
 
    \param number_of_ev number of eigen functions to read
  */
  void read_SLEPC_solution(unsigned int number_of_ev);


  //! Transform the eigenstates with S^-1/2, if needed
  void transform_eigenstates(void);

 
  //!creates constraints
  // void make_constraints(void);

 
  //!number of nodes used in the model
  unsigned int number_of_nodes;


  //! compares eigenstate energy for electrons 
  static bool compare_eigen_energy_electrons(const double& state1, const double& state2);

  //! compares eigenstate energy for holes 
  static bool compare_eigen_energy_holes(const double& state1, const double& state2);
  

  //! solutions of the eigenvalue problem
  //std::vector<eigen_propblem_solution> solution;
  

  //!list of periodic nodes
  std :: vector< std :: vector <const Node*> >  nodes_periodic; //dim node list's: each contains list of nodes that periodic b.c
                                                                //must be applied to
  

  void calculate_Hamiltonian_and_S(void);


  //! calculates the norm of the eigenstate \f$ \sqrt {| \langle \psi|\psi \rangle |} \f$
  /*!
    \param state_number number of the eigenstate
  */
  double eigenstate_norm(unsigned int state_number);


  //!Calculates Fermi Dirac probability
  /*!
    For electrons:  \f$ p = \frac{1}{1 + \exp (\frac{E - \mu}{kT})}     \f$
    For holes:      \f$ p = 1 - \frac{1}{1 + \exp (\frac{E - \mu}{kT})} \f$
    \param Energy   state energy [eV]
    \param Fermi_energy  Fermi energy [eV]
    \param temperature   temperature [K] 
   */
  double Fermi_statistics_probability(double Energy, double Fermi_energy, double temperature) const;


  //!Calculate number of bands in the Hamiltonian 
  short calculate_number_of_bands(void) const;



  //!calculated density
  std::map<const Elem*, double> _density;
 
 

  //!k-vector in atomic units
  double k_vector[3];



  //!put spectrum shift energy to be almost equal to the 1st eigenvalue
  double get_new_spectrum_shift(void);


  //!return electric potential (wrapper)
  void  get_electric_potential(const Elem* elem, const std::vector<Point>& q_point, std::vector<double> electric_potential) const;

  //!returns band edge (without potential)
  double get_band_edge( const Elem* elem) const;

  //!return electric potential from drift-diffusion
  double get_electric_potential(const Elem* elem, const Point&  qp) const;

  //!return electro chemical potential from drift-diffusion
  double get_electro_chem_potential(const Elem* elem) const; 


  //!point for bulk dispersion
  Point _bulk_point;
  


  //! bulk eigenstates
  void solve_bulk(void);

 protected:


  virtual void get_solution_secure(const Elem* elem,
      std::map<ID, std::vector<double> >& values,
      const std::vector<Point>& points);

  virtual void get_solution_secure(
      std::map<ID, std::vector<double> >& values);


  virtual void do_setup_solution_variables(void);


  virtual void 	do_init(void);

  virtual void 	do_solve (void);

  virtual void 	parse_options (void);


  //! We override this to write in our own format
  virtual void plot_globaldata(void);

 
  virtual void do_assemble(const ModelOptions& options);


};
//-------------------------------------------------------------------
inline double EnvelopFunctionApprox::Fermi_statistics_probability(double Energy, double Fermi_energy,
                                                                  double Temperature) const
{
  

  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double el_fermi;

  //if (exp_arg > 20)
  //  el_fermi = 0.0;
  //else
    el_fermi = 1.0/(  1 +  std::exp(exp_arg)  );


 
  if (opt.particle == "el")
    return(el_fermi);
  else
    return(1.0 - el_fermi);


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

inline  unsigned int EnvelopFunctionApprox::get_number_of_bands(void) const
{
  return (opt.number_of_bands);
}

//--------------------------------------------------------------------
inline const std::map<short, short>& EnvelopFunctionApprox::get_kp_bands(void) const
{

  return (opt.kp_bands);

}

//--------------------------------------------------------------------
inline double EnvelopFunctionApprox::get_particle_charge(void) const
{
  if (opt.particle == "el")
    return -1.0;
  else
    return +1.0;
} 

inline const Elem*  EnvelopFunctionApprox::return_bulk_element(void) const
{
  return _bulk_mat_element;
}

#endif
