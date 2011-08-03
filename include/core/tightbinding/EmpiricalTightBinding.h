#ifndef _EMPIRICALTIGHTBINDING_H_
#define _EMPIRICALTIGHTBINIDNG_H_

#include "TightBinding.h"

class ETB : public TightBinding
{

 public:

  enum Solutions
  {
    MeshStates,  //Eigenstate Magnitude projected on mesh
    ElQuantumDensity,  //Electron charge density
    HlQuantumDensity  //Hole charge density
  };

  typedef enum {JVXL=1, CUBE=2} OutputFormat;

  //! A class for Dftb options
  class UptOptions
  {
  public:

    UptOptions(void);
    ~UptOptions(void);

    int verbose;
    int max_TB_order;
    bool harrison_flag;
    bool relat_flag;
    bool potential_flag;
    bool band_shift_flag;
    bool opt_flag;
    bool check_bondmap;
    int poldir;
    std::string potential_sim;
    std::string strain_sim;
    std::string etb_dataset;
    char* database_path;
    char* work_path;
    char* out_path;
    char* load_path;
    char* upt_filename;
    char* gen_outfile;
    char* sparse_fmt;
    OutputFormat out_format;
    double c_axis[3];
    double k_point[3];
    double hl_chem_pot;
    double el_chem_pot;
    double temperature;
    double projection_length;
    double dg_scale;
    double dg_onsite;
    double grid_step;
    bool hybrid_passivation;
    bool d_states_correction;
  };

  class UptSolverOptions
  {
  public:
    UptSolverOptions(void);
    ~UptSolverOptions(void);

    std::string solver;
    int n_vb;
    int n_cb;
    int start_vb;
    int start_cb;
    int min_iter;
    int long_iter;
    int max_iter;
    int m0;
    double guess_vb;
    double guess_cb;
    double fast_tol;
    double long_tol;
    double ort_tol;
    double e_min;
    double e_max;
    int twice_vb;
    int twice_cb;
    bool read_states;
  };

  //! Constructor
  ETB(const ModelOptions& options);

  //! Destructor
  ~ETB(void);

  //! Create object
  static ETB* create(const ModelOptions& options);

  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);

  //! initialize or reinitialize the library container's with structure data
  void reinit(void);

  //! assemble the matrix again w/o init overheads
  void reassemble(void);

  //void set_kpoint(void);

  void set_num_states(int num_vb, int num_cb);

  //! computes the solution
  void solve_for_particle(const std::string& particle);

  //! Computes the Fermi level averaged according to the state density
  double calculate_fermi_averaged(unsigned int i);

  //! tells if the ETB calculation contains SO coupling
  bool is_relativistic(void);

  //! compute atomic charges
  void compute_atomic_charges(const std::string& particle, std::vector<double>& qmat);

  //! compute state density for a single state
  void compute_eigenvector_mag(unsigned int, std::vector<double>&);

  //! Provide solution values
  virtual void
    get_solution_secure(const Elem* elem, std::map<ID, std::vector<double>>& values,
        const std::vector<Point>& p);


 protected:

//  void build_statedens(std::vector<double>& values, const Point& r);

  double build_rho3d(const std::vector<double>& tb_density, const Point& r);

  double build_rho2d(const std::vector<double>& tb_density, const Point& r);

  double build_average_rho1d(const std::vector<double>& tb_density, const Elem* elem);

  virtual void do_init(void);

  virtual void do_solve (void);

  virtual void plot_globaldata (void);

  virtual void parse_options(void);

  virtual void do_assemble(const ModelOptions& options);

  //! Setup the available variables
  virtual void do_setup_solution_variables(void);


  //Mesh dimension (used many times by charge projection function)
  int _dim;

  //! computes the matrix element for the optical matrix
  //! the P matrix comes out in   eV * Ang  (setting m=1,hbar=1)
  virtual std::complex<double> calculate_matrix_element(const std::string& i_particle,
							unsigned int i,
							const std::string& j_particle,
							unsigned int j);

  /*! Note: for the moment calculate_matrix_element relays on the fact that the first
     *  n_vb states are valence states, then there are all the electron states.
     *  This should change if we want to put a self-consistent density calculation
     *  that will dynamically append new states to the _solution vector.
     */


 private:

  //! Get options suited for DFTB+ tight binding builder and solver
  void get_upt_options(void);

  //! Function for building options from input
  void build_input_options(void);

  //! Print all _dftb_options for debugging purposes
  void print_upt_options(void);

  void read_kpoints(void);

  //! Add potential shifts
  void add_pot_shifts(void);

  //! Add band shifts
  void add_band_shifts(void);

  void write_shifts(void);

  //! project strain on atom (Macrostrain)
  void project_atom_strain(void);

  //! subroutine used to read band-edges from database
  void get_band_edges(void);

  //! get the band extrema
  void get_band_extrema(double& cb_min, double& vb_max);

  //! get c-axis orientation
  void get_c_axis(void);

  //! Structure containing options for DFTB+ tight binding builder
  UptOptions _upt_options;

  //! Structure containing options for DFTB+ tight binding solver
  UptSolverOptions _upt_solver_options;

  //! Uptight instance associated to the simulation
  UptWrapper* inst;

  //! flag to decide whether init the structure (e.g. if strained)
  int _init;

  //! flag to decide whether assemble the matrix again on not
  int _assemble;

  //! flag to decide whether to read a structure from file
  std::string _upg_filename;

  //! vector to hold number of orbital per ion
  std::vector<int> _ion_num_orbitals;

  //! State density on atom, first index is the eigenvector label
  std::map<unsigned int, std::vector<double>> _eigenvector_mag;
 
  //! Electron charge density on atoms
  std::vector<double> _el_atomic_charges;
  
  //!Hole charge density on atoms
  std::vector<double> _hl_atomic_charges;

  //!Number of atoms (without including hydrogens)
  unsigned int _N_without_H;
  
  std::map<ID, double > _map_ID_Evb;
  std::map<ID, double > _map_ID_Ecb;  

  //!Vector for atom-projected band shifts
  std::vector<double> _band_shift;
  
  double _vb_shift;

    
  //! Size of the solution (number of states)
  unsigned int _solution_size;
  

};

// Inline members definition
//----------------------------

inline
ETB* ETB::create(const ModelOptions& options)
{
  return new ETB(options);
}

inline
void ETB::reassemble()
{
  _assemble = 1;
}

inline
void ETB::set_num_states(int num_vb, int num_cb)
{
  _upt_solver_options.n_vb = num_vb;
  _upt_solver_options.n_cb = num_cb;
}

inline
bool ETB::is_relativistic(void)
{
  return _upt_options.relat_flag;

}




#endif
