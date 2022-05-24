// $Id: Negf.h 2964 2011-10-10 20:34:57Z fpalomba $

#ifndef _NEGF_H_
#define _NEGF_H_

#include "SimulationInterface.h"
#include "EigenvalueProblem.h"
#include "AtomisticStructure.h"
#include "ModelOptions.h"
#include "QuantumContact.h"
#include "TiberLinearSystem.h"
#include "Boundary.h"
#include "KspaceIntegration.h"

#include <string>

class NegfWrapper;

/*!
 *
 * \brief This is an example implementation of the Poisson equation to
 *        help module development.
 *
 * Illustrates the basic usage of the SimulationInterface API.
 */
class TBDLLOCAL Negf : public SimulationInterface
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
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Negf(void);

    //! We need a public static creator function
    static Negf* create(const ModelOptions& options);

    void reorder(void);

    //void set_stuff_for_etb(void);

    //! call-back method that KspaceIntegration invokes
    void calculate_for_k_point(const Point& kpoint, DofField& spectrum, double& estimator);

  protected:

    //! The initialization
    virtual void do_init(void);

    //! Re-initialization for sweep, etc...
    virtual void do_reinit(void);

    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Solve the MyPoisson equation
    virtual void do_solve(void);

    void print_Lib(void);

    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
                                                  const Material* mat) const;


    virtual void do_setup_solution_variables(void);

    //! Used to plot global data such as current
    virtual void plot_globaldata (void);

    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);

    void init_hamil(void);

    void setup_hamil(void);

    void setup_negf(void);
    
    void finalize(void);

    void compute_current(void);

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

	double deltaE;

        bool writeLDOS;
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

    double get_band_edge(const std::string& band) const;
    double get_band_edge(SimulationInterface* model, const std::string& band, const Elem* elem) const;

    void get_mu_and_bands(std::vector<double>& Ec, std::vector<double>& Ev,
                          std::vector<double>& muN, std::vector<double>& muP);

    std::vector<double> get_ordered_solution(SimulationInterface* model, const std::string& var);

    const Boundary* get_boundary(const QuantumContact* qc);

    //void apply_dirichlet_bc(void);
    
    void init_efa_hamil(void);

    void init_etb_hamil(void);

    //void setup_sb_hamil(void); 
   
    void setup_efa_hamil(void);
    
    void setup_etb_hamil(void);

    bool is_generalized(void);

    void plot_LDOS(const std::vector<double>& energies,
        const std::vector<std::vector<double>>& ldos,
        const std::string& name_suffix);

    void occupy_LDOS(const std::vector<double>& ldos);

    void transfer_density(const std::vector<double>& density,
        const std::string& particle);

    //! Print energy resolved data to file
    void print_energy_resolved(const std::string& file,
        const std::vector<double>& energy,
        const std::vector<std::vector<double>>& data,
        const std::string& header = "") const;

    Device* _device;

    SimulationEnvironment* _env;

    std::map<const Boundary*, QuantumContact*> _qc_boundaries;

    std::map<const QuantumContact*, const Boundary*> _bd_map;
    
    std::map<const Boundary*, unsigned int> _bd_num;

    //! Map quantumContact ID to pointers
    std::map<ID, QuantumContact*> _quantum_contacts;

    Negf(const ModelOptions& option);

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
};




#endif
