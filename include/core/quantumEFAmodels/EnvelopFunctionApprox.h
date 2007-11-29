#ifndef _ENVELOPFUNCTIONAPPROX_H_ 
#define _ENVELOPFUNCTIONAPPROX_H_
//! A class that constructs Hamiltonian and S-matrix 


#include "Macrostrain.h"
#include "SimulationInterface.h"
#include "EigenvalueProblem.h"

class EnvelopFunctionApprox  : public EigenvalueProblem
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

    //double  length_scale;   //!< mesh length scale [Bohr radius]


    //bool periodicity[3];    //!< periodic boundary conditions

   

    std::string output_type; //!< output type 

   

    double spectrum_shift;    //!< shift of spectrum ised in matrix assembly[eV]

    bool  consider_strain;    //!< apply strain effect to the EFA Hamiltonian;

    bool  consider_potential; //!< apply strain effect to the EFA Hamiltonian;

    bool estimate_spectrum_shift; //!< calculate spectrum shift from band edges;
  


    bool convergent_density;//!< if true, the number of eigenstates will be increased to reach the tolerance

 
    unsigned int initial_eigenstates_number; //!< initial number of eigenstates that is used in an iterative calculation of the density


    double relative_density_tolerance; //!< stops itarations if \f$ \rho_i / \rho_{i+1} < \varepsilon    \f$, where \f$ \rho \f$ is the                                              total density 
 
    double eigen_number_increase_factor; //!< to increase number of eigenstates for the next iteration 


    bool log_output; //!< to do a lot of output on screen



    JobKind job; //!< a job to do


    double Temperature;//!<Temperature for density calculation


    bool local_occupation; //!<If a local occupation is considered 

   

  };



  struct eigen_energy
  {
    double energy; //!< eigen energy [eV]
    unsigned int global_number; //< eigen vector
  };



  //! data structure that contains eigenvalue and eigenvector
  struct eigen_propblem_solution
  {
    double eigen_energy; //!< eigen energy [eV]
    std::vector< std::complex<double>  > eigen_vector; //< eigen vector
    double Fermi_energy; //< electro-chemical potential [eV] \f$ \langle \psi |\mu({\bf r} | \psi \rangle  \f$
    double Temperature; //!< averaged temperature for the state [K] 
  };


 


 


  //!constructor
  EnvelopFunctionApprox(void);


  //!destructor
  ~EnvelopFunctionApprox();


  //!computes Hamiltonian and S matrix
  virtual void calculate_Hamiltonian_and_S(void); 



  

 

  //!writes on disk the eigenfunction
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_eigen_function(unsigned int state_number,  const std::string& filename);


  //!writes on disk the probability function \f$  \sum_i |\psi({\bf r})|^2   \f$
  /*!
    \param state_number eigenstate number
    \param filename name of file
  */
  void output_probability_function(unsigned int state_number,  const std::string& filename);


  

  //! returns a reference to solutions
  const std::vector<eigen_propblem_solution>& get_solution() const;
 
 
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
  
 
  //! claculate total density
  /*!
    \f$ \rho = \sum_i F_{\rm{Fermi}}(E_i) \f$
    
  */
  double get_integrated_probability();



  //!obtain convergent density
  /*!
    
   
  */
  void  calculate_convergent_density( void );
  


   //! calculate cell density  (in atomic units) for a single \f$ {\bf k}_{\|}\f$ vector.  
  /*!
    The nodal density reads: \f$ \rho({\bf r}) = \sum_i   |\psi_i({\bf r})|^2 F_{fermi}(E_i) \f$ 
    The cell density reads:  \f$ \rho = \frac{1}{\Omega_0} \int_{\Omega_0}  \sum_i   |\psi_i({\bf r})|^2 \, dV  F_{fermi}(E_i), \f$
    where \f$ \Omega_0 \f$ is the element volume.

   
  

  */
  void  calculate_density(void);


  //! sets opt.initial_eigestates_number
  void set_initial_eigenstates_number(unsigned int n);



  //! returns number of active cells
  unsigned int get_number_of_active_cells();


  static  EnvelopFunctionApprox* create();


  virtual PhysicalModel*
    create_physical_model(const ModelOptions& options) const
    throw (ModelErrorException);
    
   
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException);


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



  //!calculates analitycal charge density for a 1D quantum structure (e.g. quantum well)
  /*!
    \f$  \rho({\bf r}) = \frac{mkT}{2 \pi \hbar^2}  |\psi({\bf r})|^2 \ln (1 + \exp (\frac{\mu - E}{kT}) )    \f$
  */
  std::map<const Elem*,  double> estimate_density1D(unsigned int state_number, double parallel_mass);


  //!calculates analitycal charge density for a 1D quantum structure (e.g. quantum wire)
  /*!
    \f$  \rho({\bf r}) = |\psi({\bf r})|^2 \frac{1}{2} \sqrt{\left( \frac{mkT}{2\pi\hbar^2} \right)}  
    F_{-1/2}\left(  \frac{\mu - E}{kT}         \right)       \f$
  */
  std::map<const Elem*, double> estimate_density2D(unsigned int state_number, double parallel_mass);


  inline double get_particle_charge(void) const; 
  

  

 private:

  //!pointer to the device object
  static  Device* _device;


 
  //!pointer to mesh of the equation systems
  //Mesh* mesh;

  options opt;

  //!pointer to meshdata of the equation systems
  MeshData*  meshdata;

  //!pointer the equation systems object
  //EquationSystems* es;

  
  std::string system_name;

  //!pointer to a drift-diffusion object that is used to get potential data 
  SimulationInterface* poisson_equation;


  //!pointer to a temperature simulation object that is used to get temperature data
  SimulationInterface* temperature_simulation;

  //! ID of temperature in heat equation
  ID temperature_ID;


  //! ID of electric potential in poisson_equation
  ID potential_ID;

  //! ID of electro chemical potential in poisson_equation for particles
  ID electro_chem_pot_ID;

  //! band edge ID in poisson_equation
  ID band_edge_ID;

  //!pointer to the macrostrain object that is used to get strain data 
  Macrostrain* strain;

  //!system that we add to the equation systems
  //LinearImplicitSystem* system;

  //!diriclet nodes vector
  // std::vector<unsigned int>  dirichlet_nodes;

 
  
  //!my Jacobian because I calculate everything in atomic units
  //double my_Jacobian; 

  
 


  //!calculates \f$ |\langle \psi|\psi \rangle|^2 \f$ of an eigenstate
  /*!
    \param state_number eigenstate number
    \param prob_data \f$ |\langle \psi|\psi \rangle|^2 \f$ of an eigenstate
  */
  void prepare_probability_function(const unsigned int state_number, std::vector<double>& prob_data);






  //!bands names
  std::vector<std::string> psi_name;


 

 


  //!passes S matrix to the eigensolver
  void copy_S_matrix_to_solver();

 


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



 
 
  



  //!creates constraints
  // void make_constraints(void);

 

 


  //!number of nodes used in the model
  unsigned int number_of_nodes;


  //! compares eigenstate energy for electrons 
  static bool compare_eigen_energy_electrons(eigen_propblem_solution state1, eigen_propblem_solution state2);

  //! compares eigenstate energy for holes 
  static bool compare_eigen_energy_holes(eigen_propblem_solution state1, eigen_propblem_solution state2);

  
  //! compares eigenstate energy for electrons 
  static bool compare_eigen_energy_electrons1(eigen_energy state1, eigen_energy state2);

  //! compares eigenstate energy for holes 
  static bool compare_eigen_energy_holes1(eigen_energy state1, eigen_energy state2);
  

  //! solutions of the eigenvalue problem
  std::vector<eigen_propblem_solution> solution;
  


 

  //! Apply periodic boundary conditions
  //void apply_periodic_bc();

  //! create list of nodes that lies at the periodic boundary
  //void make_nodes_periodic();


  //!list of periodic nodes
  std :: vector< std :: vector <const Node*> >  nodes_periodic; //dim node list's: each contains list of nodes that periodic b.c
                                                                //must be applied to
  


  //! calculates the norm of the eigenstate \f$ \sqrt {| \langle \psi|\psi \rangle |} \f$
  /*!
    \param state_number number of the eigenstate
  */
  double eigenstate_norm(unsigned int state_number);

 
  
  //! calculate density \f$ | \psi_i (r) |^2 \f$
  /*!
    \param i number of the eigenstate
   */
  std::vector<double> calculate_prob_function(unsigned int i);


  //! calculate density  \f$ \frac{1}{\Sigma_0} \int_{\Sigma_0} | \psi_i (r) |^2 /,dV \f$
  /*!
    \param i number of the eigenstate
  */
  std::vector<double> calculate_cell_prob_function(unsigned int i);


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


  //!in this class  outputs \f$ |<\psi|\psi>|^2 \f$ for each eigenstate
  virtual void build_nodal_results(const std::set<std::string>& variables,
				   std::vector<double>& results, std::vector<std::string>& legend);


  virtual void 	build_elemental_results (const std::set< std::string > &variables, 
					 std::vector< double > &results, std::vector< std::string > &legend);


  //!in this class it outputs eigen energies
  virtual void 	build_integrated_quantities (const std::set< std::string > &names, std::vector< double > &values);

  //!in this class it builds descriotion for eigen energies
  virtual void 	build_integrated_quantities_description (const std::set< std::string > &names,
							 std::vector< std::string > &legend, 
							 std::vector< std::string > &description);
 



  //!in this class outputs eigenvalues 
  void get_integrated_quantities(const std::set<std::string>& names,
				 std::vector<double>& values);


 


  virtual void 	do_init(void);

  virtual void 	do_solve (void);

  virtual void 	parse_options (void);


 



};
//-------------------------------------------------------------------
inline double EnvelopFunctionApprox::Fermi_statistics_probability(double Energy, double Fermi_energy,
                                                                  double Temperature) const
{
  

  double T_EV = Temperature * Constants::k_Boltzmann;
  double exp_arg =  (Energy - Fermi_energy)/T_EV;
  
  double el_fermi;

  if (exp_arg > 20) 
    el_fermi = 0.0;
  else
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

inline EnvelopFunctionApprox*  EnvelopFunctionApprox::create()
{
  return (new EnvelopFunctionApprox );
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



#endif
