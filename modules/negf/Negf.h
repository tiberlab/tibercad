/*  
 * This file is part of the tiberCAD module negf.
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
 * \file Negf.h
 * \brief tiberCAD negf module header.
 *
 * \note This file is part of module negf.
 */


#ifndef TC_NEGF_H
#define TC_NEGF_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/physics/schroedinger/EigenvalueProblem.h"
#include "tibercad/atomistic/AtomisticStructure.h"
#include "tibercad/base/ModelOptions.h"
#include "tibercad/geom/QuantumContact.h"
#include "tibercad/solver/TiberLinearSystem.h"
#include "tibercad/geom/Boundary.h"
#include "tibercad/kintegration/KspaceIntegration.h"
#include "tibercad/kintegration/Kspace.h"

#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"

#include <string>

class NegfWrapper;

/*!
 *
 * \brief This is an example implementation of the Poisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class TC_DLLOCAL Negf : public SimulationInterface
{

  public:

    enum Solutions
    {
       ReorderPotential,  // Laplace equation solution from reordering routines
       elDensity,  // Electron QuantumDensity from Negf
       hlDensity,  // Hole QuantumDensity from Negf
       eCurrentDensity,    //
       hCurrentDensity,    // Contact Currents
       LDOS,               // LDOS and occupation
       ContactCurrent = 100
    };

    enum kIntegrationType
    {
       INTDENSITYEL,
       INTDENSITYHL, 
       INTCURRENT
    };

    //! Destructor
    virtual ~Negf(void);

    void reorder(void);

    //! call-back method that KspaceIntegration invokes
    void calculate_for_k_point(const Point& kpoint, DofField& spectrum, double& estimator);

  protected:

    //! Constructor
    explicit Negf(const ModelOptions& option);

    //! The initialization
    virtual void do_init(void) override;

    //! Re-initialization for sweep, etc...
    virtual void do_reinit(void) override;

    //! Parse the options from the input file
    void parse_options(void);


    //! Solve the MyPoisson equation
    virtual void do_solve(void) override;


    //! Reimplement the MPI communicators setup (TEMPORARILY COMMENTED UNTIL API IS FINISHED)
    virtual void setup_mpi_comm(void) override;


    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
                                                  const Material* mat) const override;


    virtual void do_setup_solution_variables(void) override;

    //! Used to plot global data such as current
    virtual void plot_globaldata (void) override;

    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p) override;

    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values) override;

    void init_hamil(void);

    void setup_hamil(void);

    void setup_negf(void);
    
    void finalize(void);

    void compute_current(void);
    
    std::vector<double> compute_layer_currents(void);

    void reorder_assemble(void);

    void get_boundary_potentials(QuantumContact* qc, double& av_V, double& av_mu_n, double& av_mu_p);

    //void calculate_density(const std::string& particle);

    void init_k_space_integration(void);
    
    void init_k_space(ModelOptions& opt);

    struct options
    {
        std::string pot_module;

        double Emin;

        double Emax;

        double Estep;

        double Estep_coarse;

        std::vector <int> Np_n;

        int n_poles;

        int n_kT;

        double Np_real;

        double DEc;

        double DEv;

        int verbosity;

        double delta;

        double dos_delta;

	      double deltaE;

        bool writeLDOS;

        bool contact_calculation;

        bool scattering;
        
        //! discretization along the transport direction needed for computation of the elph coupling
        double deltaz;
        
        //! area of the cell for computing inelastic coupling
        double cell_area;

        //! volume of the cell for computing inelastic coupling
        double volume;

        // bool dump_matrices;
    };




  private:

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
      MyAssembly(Negf* obj) : _obj(obj) {};

      void assemble() override
      {
        _obj->reorder_assemble();
      }

      private:
      Negf* _obj;

    };

    MyAssembly _reorder_assembly;

    enum InteractionModels
    {
      DUMMY,
      DEPHDIAGONAL,
      DEPHATOMBLOCK,
      DEPHOVERLAP,
      POLAROPTICAL,
      NONPOLAROPTICAL,
      ACOUSTICINEL,
      PHOTON,
    };

    //! The structure containing elastic and inelastic scattering parameters. 
    //! TODO: Some of these fields could be unnecessary
    struct Interaction
    {
      // Specify which model in input
      int model = DUMMY;

      // Coupling strength
      double coupling = 0.0;

      // Iterations for self-consistent Born approximation
      int scba_niter = 0;

      // Tolerance for self-consistent Born approximation
      double scba_tol = 1e-7;

      // List of orbitals per atom
      std::vector<int> orbsperatm = std::vector<int>(1, 0);

      // phonon frequency
      double wq = 0.0;

      // dielectric constant in the high freq limit
      double eps_inf = 1.0;

      // dielectric constant in the low freq
      double eps_r = 1.0;

      // screening paramter
      double q0 = 1.0;

      // deformation potential
      double D0 = 0.0;

      // whether tri-diagonal blocks are computed
      bool tTridiagonal = true;

      // the intensity of incident light
      double intensity = 0.0;

      // the energy of the incident photon
      double Ephot = 0.0;

      // the refraction index
      double nr = 1.0;

      // the direction of the polarization vector (1-indexed)
      int poldir = 3;

      int deb_id = 0;
    };

    double get_band_edge(const std::string& band) const;
    double get_band_edge(SimulationInterface* model, const std::string& band, const Elem* elem) const;

    void get_mu_and_bands(std::vector<double>& Ec, std::vector<double>& Ev,
                          std::vector<double>& muN, std::vector<double>& muP);

    std::vector<double> get_ordered_solution(SimulationInterface* model, const std::string& var);

    const Boundary* get_boundary(const QuantumContact* qc);

    //void apply_dirichlet_bc(void);
    
    void init_efa_hamil(void);

    void init_etb_hamil(void);

    void reinit_etb_module(void);

    //void setup_sb_hamil(void); 
   
    void setup_efa_hamil(void);
    
    void setup_etb_hamil(void);

    bool is_generalized(void);

    void plot_LDOS(const std::vector<double>& energies,
        const std::vector<std::vector<double>>& ldos,
        const std::string& name_suffix);

    void occupy_LDOS(const std::vector<double>& ldos);

    void create_python_script(void);

    void transfer_density(const std::vector<double>& density,
        const std::string& particle);

    void transfer_density_efa(const std::vector<double>& density,
        const std::string& particle);

    void transfer_density_etb(const std::vector<double>& density,
        const std::string& particle);

    double project_density(const Elem* elem, const Point& point, 
        const std::vector<double>& atomic_charges, double cutoff);

    //! Print energy resolved data to file
    void print_energy_resolved(const std::string& file,
        const std::vector<double>& energy,
        const std::vector<std::vector<double>>& data,
        const std::string& header = "") const;

    void get_equivalent_points(Kspace* kspace, const std::vector<DofField> kpoints, 
         std::vector<DofField>& equiv_points, std::vector<int>& n_equiv);

    //! Wrapper for setting the kpoints that depends on the integrated quantity
    void set_kpoints(std::string solution);

    //! Passes kpoints to libNEGF, after calculating global and local indices, either of the full BZ 
    //! or the reduced one (with equivalent kpoints)
    void set_kpoints(KspaceIntegration* k_int, bool reduced_BZ);

    //! Used within set_kpoints and get_polarization_matrix 
    void get_distributed_kpoints(KspaceIntegration* k_int, std::vector<double>& kweights, std::vector<int>& local_k_indices, 
                                            std::vector<std::vector<double>>& global_abs_kpoints, bool reduced_BZ);
    
    //! Passes the Hamiltonians to libNEGF for each kpoint
    void set_hamiltonians(void);

    //! Computes (R_a - R_b) * H_ab(k) for each k (in CSR format)
    // void get_polarization_matrices(std::vector<std::vector<std::vector<std::complex<double>>>>& P, int poldir);
    void get_polarization_matrices(std::vector<std::vector<int>>& IP, std::vector<std::vector<int>>& JP, 
                                   std::vector<std::vector<std::complex<double>>>& P, int poldir);

    //! Pass information about atomistic structure such as: lattice vectors, coordinates of atoms, matrix_indices.
    void init_basis(void);

    std::vector<DofField> transform_to_fractional_coordinates(Kspace* kspace, const std::vector<DofField>& kpoints);

    void get_coordinates(std::vector<DofField>& coordinates);

    std::vector<DofField> get_lattice_vectors(std::vector<double>& r1, std::vector<double>& r2);
    
    std::vector<DofField> get_lattice_vectors(std::vector<double>& r1);

    //! Parse the options for interactions from the input file.
    void parse_scattering_options(void);

    int get_scattering_model(std::string model);

    void setup_interactions(void);

    void print_interactions(void);

    void  compute_area_and_volume(void);

    template<typename T> void compute_area_and_volume(const T& a, const T& b, const T& c);

    Device* _device;

    SimulationEnvironment* _env;

    std::map<const Boundary*, QuantumContact*> _qc_boundaries;

    std::map<const QuantumContact*, const Boundary*> _bd_map;
    
    std::map<const Boundary*, unsigned int> _bd_num;

    //! Map quantumContact ID to pointers
    std::map<ID, QuantumContact*> _quantum_contacts;

    //! This system is used for calculating DOF ordering
    TiberLinearSystem* _sys;

    //! This system is used to work with DOF indices
    TiberLinearSystem* _sys_H;

    //! This system is used to handle density results
    TiberLinearSystem* _qdens_sys;

    //! permutation vector (for DOFs)
    std::vector<unsigned int> _perm;

    //! inverse permutation vector (for DOFs)
    std::vector<unsigned int> _inv_perm;

    //! The ends of the PL blocks
    std::vector<int> _end_blocks;

    //! The ends of the contact blocks
    std::vector<int> _contend;

    //! The start of the surface blocks
    std::vector<int> _surfstart;

    //! The ends of the surface blocks
    std::vector<int> _surfend;

    //! The indices of the contact blocks
    std::vector<int> _cblk;

    unsigned int _device_n_dofs;

    std::vector<unsigned int> _qc_n_dofs;

    NegfWrapper* _libnegf;

    options opt;


    KspaceIntegration* _k_int_density;

    KspaceIntegration* _k_int_current;

    std::vector<Interaction> _interactions;


    //!k-dependent index for setting the Hamiltonian in the HS container
    int _iK;

    //! local number of Hamiltonians to be passed to libNEGF
    int _n_Hk;

    //! kpoints broadcasted to all processes, in absolute units (inverse length)
    std::vector<DofField> _global_abs_kpoints;
    
    //! local indices that map the local kpoints to the global ones 
    std::vector<int> _local_k_indices;

    //! the lattice vectors to be passed to init_basis
    std::vector<DofField> _lattice_vectors;

    // internal status for k_space_integration
    int _which_integration;

    libMesh::VectorValue<double> _k_vec;

    DofField current;

    std::map<const QuantumContact*, double> _contact_potential;

    std::map<const QuantumContact*, double> _contact_current;

    double mumin, mumax;

    std::string _hamil_type;

    //! The provider of the hamiltonian
    EigenvalueProblem* _ext_module;

    AtomisticStructure* _atom_structure;

    // The MPI cartesian communicators used in libNEGF
    libMesh::Parallel::Communicator _cart_comm;
    libMesh::Parallel::Communicator _k_comm;
    libMesh::Parallel::Communicator _en_comm;
};




#endif
