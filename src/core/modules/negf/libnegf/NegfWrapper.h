#ifndef _NEGFWRAPPER_H_
#define _NEGFWRAPPER_H_

#include "libnegf.hpp"
#include "lnParams.h"

#include <complex>
#include <vector>
#include <map>

#include <mpi.h>
//-----------------------------------------------------------------------
typedef std::complex<double> Complex; 

//!Wrapper class for callings to libNEGF library (libnegf.so)
class NegfWrapper
{

public:

  typedef struct lnparams Parameters;

  //!Constructor
  /*!Assign an handler to libNEGF instance, transparent to programmer
   *
   */
  NegfWrapper(void);

  //! Destructor  
  ~NegfWrapper(void);

  //!Static method to create a Negf wrapper instance
  static NegfWrapper* create(void);

  //!Initialize the Negf instance
  void init(void);

  //Initialize MPI communicator
  int set_mpi_comm(MPI_Comm comm);

  //Initialize MPI k-E cartesian communicators
  void mpi_cart_init(MPI_Comm inComm, int nGroups, MPI_Comm& cartComm, MPI_Comm& kComm, MPI_Comm& enComm);

  //! Destroy and reinitialize
  void force_reinit(void);

  //!read H and S hamiltonians 
  int read_HS(std::string real_path, std::string imag_path, std::string H_or_S);

  //!read 'negf.in'
  int read_input(void);

  //! Get the parameters
  void get_parameters(NegfWrapper::Parameters& params);

  //! Set the parameters
  void set_parameters(const NegfWrapper::Parameters& params);

  //!Compute current
  //! unitsOfH is a string like "H", "eV", "cm^-1"
  //! unitsOfJ is a string like "A"
  double current(std::string unitsOfH, std::string unitsOfJ);

  //! Compute the currents flowing through each layer
  //! unitsOfH is a string like "H", "eV", "cm^-1"
  //! unitsOfJ is a string like "A"
  void layer_current(std::vector<double>& layer_current, std::string unitsOfH, std::string unitsOfJ);

  //! Get the energy points
  /*!
   * Only real part of energies is given back. These can be used e.g.
   * to print transmission or LDOS
   */
  void get_energies(std::vector<double>& energies);

  //! Return transmission
  void get_transmission(std::vector<std::vector<double>>& transmission);

  //! Return transmission
  void get_energy_current(std::vector<std::vector<double>>& current);

  //!Compute charge density
  void density(std::vector<double>& density, std::string particle, bool contact_calc);

  //!Compute charge density, using quasi-equilibrium approximation
  void quasi_equilibrium_density(std::vector<double>& density, std::string particle,
                                 std::vector<double>& Ec, std::vector<double>& Ev,
                                 std::vector<double>& muN, std::vector<double>& muP);

  //! Return the LDOS in matrix form
  void get_ldos(std::vector<std::vector<double>>& ldos);

  //!Set SC iteration
  void set_iteration(int iter);

  //!Set scratch folder path
  void set_scratch_path(std::string path);

  //!Set scratch folder path
  void set_output_path(std::string path);

  //!Clean the Negf instance variable space
  void clean_libnegf(void);

  //!Get libNEGF instance handler
  inline const int* get_handler(void){ return _handler; };

  //! Set verbosity level for the library screen output
  void set_verbose(int verbose_lev);

  //! used to set an integer value to label tunneling files
  void set_kpoint(int kpoint);

  //! Set global and local kpoints, along with equivalent kpoints for extending the irreducible wedge
  void set_kpoints(std::vector<std::vector<double>> kpoints, std::vector<double> kweights,  
       std::vector<int> local_k_indices, std::vector<std::vector<double>> equivalent_points, 
       std::vector<int> n_equivalent, bool reduced_BZ);

  //! Pass atom coordinates and map from atom indices to matrix indices, as well as lattice vectors
  void init_basis(std::vector<std::vector<double>> coordinates, std::vector<int> matrix_indices, 
       std::vector<std::vector<double>> lattice_vectors, int transport_direction);

  //! used to set the reference contact 0= min mu; 1= max mu
  void set_reference(int minmax);

  //! Set what to do with outer parts 0=none; 1=upper block; 2=full
  void device_contact_dm(int outer);

  //! Initialize structure to hold LDOS
  /*!
   * The two arrays specify the start and end indices for the
   * LDOS projection.
   */
  void init_ldos(const std::vector<int>& start,
                 const std::vector<int>& end);

  //! Initialize structure to hold LDOS
  /*!
   * In this version the number of DOFs has to be passed, so that
   * the entire LDOS matrix will be passed back.
   */
  void init_ldos(unsigned int nldos);

  //! Initialize structure to hold LDOS
  /*!
   * Here the DOF indices for a certain LDOS projection are given.
   * init_ldos(int) needs to be called first
   */
  void set_ldos_indices(unsigned int dos_index,
                        const std::vector<int>& indices);

  //! Initialized the contacts
  void init_contacts(int n_cont);

  //! Initialize structure
  void init_structure(int ncont, const std::vector<int>& surfstart,
    const std::vector<int>& surfend, const std::vector<int>& contend,
    int npl, const std::vector<int>& plend, const std::vector<int>& cblks);

  //! Find indices contact blocks
  std::vector<int> contact_blocks(int ncont, const std::vector<int>& surfstart,
      const std::vector<int>& surfend, const std::vector<int>& contend, int npl,
      const std::vector<int>& plend);

  //! Get library version
  void get_version(void);

  //! write partition info 
  void partition_info(void);
  
  //! set H via memcopy to libnegf container
  void set_H_csr(int nrow, std::vector<Complex >& A, std::vector<int>& JA,
                           std::vector<int>& IA, int iK);

  void create_HS_container(int n_Hk);

  //! set S via memcopy to libnegf container
  void set_S_csr(int nrow, std::vector<Complex >& A, std::vector<int>& JA,
                           std::vector<int>& IA, int iK);
  //! Set S as and Identity matrix
  void set_S_id(int nrow, int iK);

  //! Print matrix
  void print_mat(void);

  //! Print memory statistics
  void print_memory_statistics(void);

  void set_elph_dephasing(std::vector<double>& coupling, int scba_niter);

  void set_elph_block_dephasing(std::vector<double>& coupling, std::vector<int>& orbsperatm, int scba_niter);

  void set_elph_s_dephasing(std::vector<double> &coupling, std::vector<int>& orbsperatm, int scba_niter);

  void set_elph_polaroptical(std::vector<double>& coupling,  double wq, double kbT, double deltaz, double eps_r,  double eps_inf, 
                             double q0, double cell_area, int scba_niter,  bool tTridiagonal);

  void set_elph_nonpolaroptical(std::vector<double>& coupling, double wq, double kbT, double deltaz, double D0, double cell_area, 
                                int scba_niter, bool tTridiagonal);

  void set_elphot(std::vector<double>& coupling, double cell_vol, double nr, double Ephot, double intensity, std::vector<std::vector<int>>& IP,
                  std::vector<std::vector<int>>& JP, std::vector<std::vector<Complex>>& P, int poldir, int scba_niter, bool tTridiagonal, int deb_id);

  void set_scba_tolerances(double elastic_tol, double inelastic_tol);

private:
  int _handler[NEGF_HSIZE];

  bool _is_initialized;

};


#endif
