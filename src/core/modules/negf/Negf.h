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
#include "libnegf/NegfWrapper.h"
#include "KspaceIntegration.h"

#include <string>

struct sortclass{
  sortclass(const std::vector<Atom>& atoms) : _atoms(atoms) {}
  ~sortclass(){};
  bool operator() (int i, int j) { return (_atoms[i].get_position()(0)<_atoms[j].get_position()(0)); }
  const std::vector<Atom>& _atoms;
};

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

    void do_assemble(const ModelOptions& opt);

    void print_ham(std::string form);

    void print_Lib(unsigned int n_vars, double Ec, double Ev);

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

    static void reorder_assemble(EquationSystems& es, const std::string& system_name);

    static void ham_assemble(EquationSystems& es, const std::string& system_name);

    static Negf* static_this;

    static bool compare(ID i, ID j);

    void do_reorder_assemble(EquationSystems& es, const std::string& system_name);

    void do_ham_assemble(EquationSystems& es, const std::string& system_name);

    bool do_compare(ID i, ID j);

    void get_boundary_potentials(QuantumContact* qc, double& av_V, double& av_mu_n, double& av_mu_p);

    void calculate_density(const std::string& particle);

    void init_k_space_integration(void);
    
    void init_k_space(ModelOptions& opt);

    struct options
    {
        std::string pot_module;

        double Emin;

        double Emax;

        double Estep;

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

        bool set_dirichlet_bc;
    };



  private:

    void activate_quantum_contacts(void);
    void deactivate_quantum_contacts(void);

    double get_band_edge(const std::string& band) const;
    double get_band_edge(SimulationInterface* model, const std::string& band, const Elem* elem) const;

    const Boundary* get_boundary(const QuantumContact* qc);

    void apply_dirichlet_bc(void);
    
    void init_efa_hamil(void);

    void init_etb_hamil(void);

    void setup_sb_hamil(void); 
   
    void setup_efa_hamil(void);
    
    void setup_etb_hamil(void);

    bool is_generalized(void);

    void plot_LDOS(const std::vector<double>& ldos,
        const std::string& mod = "LDOS");
    void occupy_LDOS(const std::vector<double>& ldos);

    void transfer_density(const std::vector<double>& density,
        const std::string& particle);

    Device* _device;

    SimulationEnvironment* _env;

    std::set<const Boundary*> _dirichlet_boundaries;

    std::map<const Boundary*, QuantumContact*> _qc_boundaries;

    std::map<const QuantumContact*, const Boundary*> _bd_map;
    
    std::map<const Boundary*, unsigned int> _bd_num;

    //! Map quantumContact ID to pointers
    std::map<ID, QuantumContact*> _quantum_contacts;

    Negf(const ModelOptions& option);

    TiberLinearSystem* _sys;
    TiberLinearSystem* _sys_H;
    TiberLinearSystem* _sys_S;
    TiberLinearSystem* _qdens_sys;

    std::vector<ID> _perm;      //permutation vector
    std::vector<ID> _inv_perm;     //inverse permutation

    std::vector<ID> _end_blocks;

    unsigned int _device_n_dofs;

    std::vector<unsigned int> _dev;

    std::vector<unsigned int> _qc_n_dofs;

    NegfWrapper* _libnegf;

    options opt;

    //!diriclet DOFS
    std::set<unsigned int>  dirichlet_dofs;

    KspaceIntegration* _k_int_density;

    KspaceIntegration* _k_int_current;

    // internal status for k_space_integration
    int _which_integration;

    VectorValue<double> _k_vec;

    DofField _eldensity;

    DofField _hldensity;

    DofField current;

    std::map<const QuantumContact*, double> _contact_potential;

    std::map<const QuantumContact*, double> _contact_current;

    double mumin, mumax;

    std::string _hamil_type;

    EigenvalueProblem* _ext_module;

    AtomisticStructure* _atom_structure;
};




#endif
