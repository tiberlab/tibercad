#ifndef _EMPIRICALTIGHTBINDING_H_
#define _EMPIRICALTIGHTBINIDNG_H_

#include "TightBinding.h"

class ETB : public TightBinding
{

 public:

//   enum Variables
//   {
//     UNKNOWN = 0,
//     EL_CH,
//     HL_CH
//   };

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
    char* upt_filename;
    char* gen_outfile;
    char* sparse_fmt;
    double* c_axis;
    //double vb_shift;
    double hl_chem_pot;
    double el_chem_pot;
    double temperature;
    double projection_length;
    double dg_scale;
    double dg_onsite;
  };

  class UptSolverOptions
  {
  public:
    UptSolverOptions(void);
    ~UptSolverOptions(void);

    std::string solver;
    int n_vb;
    int n_cb;
    int min_iter;
    int long_iter;
    int max_iter;
    double guess_vb;
    double guess_cb;
    double fast_tol;
    double long_tol;
    double ort_tol;

  };

  //! Constructor
  ETB(void);

  //! Destructor
  ~ETB(void);

  //! Create object
  static ETB* create();

  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);

  //! initialize or reinitialize the library container's with structure data
  void reinit(void);

  void reassemble(void);

  void set_num_states(int num_vb, int num_cb);

  //! computes the solution
  void solve_for_particle(const std::string& particle);

  //! Computes the Fermi level averaged according to the state density
  double calculate_fermi_averaged(unsigned int i);

  //! tells if the ETB calculation contains SO coupling
  bool is_relativistic(void);

  //! compute atomic charges
  void compute_atomic_charges(const std::string& particle, std::vector<double>& qmat);

  virtual void get_solution_secure(const Elem* elem,
      const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);

  virtual void
  get_solution_secure(const Elem* elem, const std::vector<Point>& p,
      const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);

  //! Order the solution in correct mode
  virtual void build_elemental_results(const std::set<std::string>& variables,
      std::vector<double>& results, std::vector<std::string>& legend);


 private:

  //! Get options suited for DFTB+ tight binding builder and solver
  void get_upt_options();

  //! Function for building options from input
  void build_input_options(void);

  //! Print all _dftb_options for debugging purposes
  void print_upt_options(void);

  void read_kpoints(void);

  //! Add potential shifts
  void add_pot_shifts(void);

  //! Add band shifts
  void add_band_shifts(void);

  //! subroutine used to read band-edges from database
  void get_band_edges(void);

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

  /*! \copydoc SimulationInterface::convert_variable_name_to_id() */
  virtual ID convert_variable_name_to_id(const std::string& variable_name) const;
  
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
    
    //! Charge density on elements (for faster scc calculation)
    std::map<const Elem*, double> _elemental_result_el;
    std::map<const Elem*, double> _elemental_result_hl;


 protected:

  double build_rho(const std::string& particle, const Point& r);

  double build_rho3d(const std::string& particle, const Point& r);

  double build_average_rho1d(const std::string& particle, const Elem* elem);

  virtual void do_init(void);

  virtual void do_solve (void);

  virtual void do_plot (void);

  virtual void parse_options(void);

  virtual void assemble(const ModelOptions& options);

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

};

// Inline members definition
//----------------------------

inline
ETB* ETB::create()
{
  return new ETB;
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
