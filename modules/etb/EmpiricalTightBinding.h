/*  
 * This file is part of the tiberCAD module empirical_tb.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file EmpiricalTightBinding.h
 * \brief tiberCAD empirical_tb module header.
 *
 * \note This file is part of module empirical_tb.
 */

#ifndef _EMPIRICALTIGHTBINDING_H_
#define _EMPIRICALTIGHTBINIDNG_H_

#include "tibercad/physics/tightbinding/TightBinding.h"
#include "tibercad/physics/StrainInterface.h"

class ETB : public TightBinding
{

 public:

  enum Solutions
  {
    MeshStates,  //Eigenstate Magnitude projected on mesh
    MeshStatesNodes,  //Eigenstate Magnitude projected on mesh nodes
    eDensity,  //Electron charge density on mesh nodes
    hDensity  //Hole charge density on mesh nodes
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
    bool pdos_flag;
    bool potential_flag;
    bool band_shift_flag;
    bool opt_flag;
    bool check_bondmap;
    int poldir;
    std::string potential_sim;
    std::string strain_sim;
    //std::string etb_dataset;
    char* database_path;
    char* work_path;
    char* out_path;
    char* load_path;
    char* write_state;
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
    bool assemble_H;
    bool compute_densities;
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
    bool read_states;
    int dynamic;
    double bitoff;
    int solver_flag;
  };

  //! Constructor
  ETB(const ModelOptions& options);

  //! Destructor
  ~ETB(void);

  //! Create object
  static ETB* create(const ModelOptions& options);

  //! assemble the matrix again w/o init overheads
  void reassemble(void);

  //void set_kpoint(void);

  void set_num_states(int num_vb, int num_cb);

  //! Computes the Fermi level averaged according to the state density
  double calculate_fermi_averaged(unsigned int i);

  //! tells if the ETB calculation contains SO coupling
  bool is_relativistic(void);

  //! compute atomic charges
  void compute_atomic_charges(const std::string& particle, std::vector<double>& qmat);

  //! compute state density for a single state
  void compute_eigenvector_mag(unsigned int, std::vector<double>&);

  //! compute projected density of states for each atom
  void compute_pdos(void);

  //! print out H in matlab format 
  void print_H(const std::string& outpath) const;

  unsigned int get_number_of_bands(void) const;

 protected:

  //! get band-edges including potential
  double get_band_edge(const std::string& band);

  double build_rho3d(const std::vector<double>& tb_density, const Elem* elem, const Point& r);

  double build_rho2d(const std::vector<double>& tb_density, const Elem* elem, const Point& r);

  double build_average_rho1d(const std::vector<double>& tb_density, const Elem* elem);

  virtual void do_init(void) override;

  //! initialize or reinitialize the library container's with structure data
  virtual void do_reinit(void) override;

  virtual void do_solve(void) override;

  virtual void do_solve_for_kpoint(const Point& k_point) override;

  virtual void plot_atomisticdata (void) override;

  virtual void parse_options(void) override;

  std::pair<unsigned int, double> read_slepc_solution(void) override;

  //! Provide solution values
  virtual void
    get_solution_secure(const Elem* elem, std::map<ID, std::vector<double>>& values,
        const std::vector<Point>& p) override;

  virtual void do_assemble(const ModelOptions& options) override;

  virtual void do_calculate_density_at_k(DofField& density) override;

  //! Setup the available variables
  virtual void do_setup_solution_variables(void) override;

  virtual void do_copy_H_to_solver(void) override;

  virtual int get_H_dim(void) const override;
  
  virtual int get_H_nnz(void) const override;
  
  virtual void get_H_csr(std::vector<libMesh::Complex>& A,
                         std::vector<int>& JA,
                         std::vector<int>& IA,
                         const std::vector<unsigned int>& perm
                               = std::vector<unsigned int>(0)) const override;

  virtual void setup_atomistic_structure(void);

  //! Implementation of projection on different BZ
  virtual void do_project_to_primitive_cell(
      const std::vector<eigen_problem_solution>& a,
      const std::vector<libMesh::Point>& kpoints,
      std::vector<std::vector<libMesh::Complex>>& weights) const override;

  //Mesh dimension (used many times by charge projection function)
  int _dim;

  //! computes the matrix element for the optical matrix
  /*!
   * the P matrix comes out in   eV * Ang  (setting m=1,hbar=1)
   */
  virtual std::complex<double> calculate_matrix_element(const std::string& i_particle,
							unsigned int i,
							const std::string& j_particle,
							unsigned int j) override;

  /* Note: for the moment calculate_matrix_element relays on the fact that the first
   *  n_vb states are valence states, then there are all the electron states.
   *  This should change if we want to put a self-consistent density calculation
   *  that will dynamically append new states to the _solution vector.
   */


 private:

  //! Overrides the method in SimulationInterface.C
  void setup_mpi_comm(void);

  //! Get options suited for DFTB+ tight binding builder and solver
  void get_upt_options(void);

  //! Function for building options from input
  void build_input_options(void);

  //! Print all _dftb_options for debugging purposes
  void print_upt_options(void);

  //! Print upg file (etb_dataset may be changed into type options)
  void print_upg(const std::string &path, bool band_offsets = false);

  //! Prepare and call uptight
  void call_uptight(void);

  void read_kpoints(void);

  //! Add potential shifts
  void add_pot_shifts(void);

  //! extract valence band shifts for pure materials from database
  void get_valence_band_shifts(const Material* mat);

  //! Add band shifts
  void add_band_shifts(void);

  void write_shifts(void);

  //! project strain on atom (Macrostrain)
  void project_atom_strain(void);

  //! Project atomic densities onto mesh
  std::pair<double, double> project_densities(
      const Elem* elem, const Point& point, double cutoff);

  //! subroutine used to read band-edges from database
  void get_bulk_edges(void);

  //! get the band extrema
  void get_band_extrema(double& cb_min, double& vb_max);

  //! set the band extrema to member values (introduced after refactoring in do_reinit())
  void set_band_extrema(void);

  //! get c-axis orientation
  void get_c_axis(void);

  //! create a dummy null H and pass it to uptight (mainly for debug)
  void create_dummy_H(void);

  //! compute the Hamiltonian dimension using uptight
  unsigned int compute_H_dim(void);

  //! get the uptight orbital ids for given names
  void get_orbital_ids(const std::vector<std::string>& names, std::set<int>& ids) const;

  //! get projection on VB and CB atomic orbitals
  /*!
   * This is temporary, I would like to implement a more generic projection
   * method with API shared with EFA
   */
  void project_on_band_orbitals(const std::vector<libMesh::Complex>& solution,
                                std::vector<double>& projections) const;

  //! Structure containing options for DFTB+ tight binding builder
  UptOptions _upt_options;

  //! Structure containing options for DFTB+ tight binding solver
  UptSolverOptions _upt_solver_options;

  //! Uptight instance associated to the simulation
  UptWrapper* inst;

  //! flag to decide whether init the structure (e.g. if strained)
  bool _init;

  //! flag to decide whether assemble the matrix again on not
  bool _assemble;

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

  // map atom pairs to shifts
  std::map<std::pair<Specie, Specie>, double> _map_pairs_Evb;

  //!Vector for atom-projected band shifts
  std::vector<double> _band_shift;
  
  double _vb_shift = 0;

  //! the CB minimum
  double _cb_min = 0;

  //! the VB maximum
  double _vb_max = 0;

  //! Size of the solution (number of states)
  unsigned int _solution_size = 0;
  
  //! strain interface
  StrainInterface _strain_int;

  //! source for the potentials
  SimulationInterface* _dd_int = nullptr;

};

// Inline members definition
//----------------------------


inline
void ETB::reassemble()
{
  _assemble = true;
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

//-------------------------------------------------------------------------

#endif
