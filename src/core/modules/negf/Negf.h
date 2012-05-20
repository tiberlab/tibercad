// $Id: Negf.h 2964 2011-10-10 20:34:57Z fpalomba $

#ifndef _NEGF_H_
#define _NEGF_H_

#include "SimulationInterface.h"
#include "ModelOptions.h"
#include "QuantumContact.h"
#include "TiberLinearSystem.h"
#include "Boundary.h"
#include "libnegf/NegfWrapper.h"
#include "KspaceIntegration.h"

#include <string>

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
       hDensity,  // Hole QuantumDensity from Negf
       eCurrentDensity,    //
       hCurrentDensity,    // Contact Currents
       ContactCurrent = 100
    };

    enum kIntegrationType
    {
       kintegrationdensity

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

    void setup_effectivemass_hamil(void);

    void compute_current(void);

    static void reorder_assemble(EquationSystems& es, const std::string& system_name);

    static void ham_assemble(EquationSystems& es, const std::string& system_name);

    static Negf* static_this;

    static bool compare(ID i, ID j);

    void do_reorder_assemble(EquationSystems& es, const std::string& system_name);

    void do_ham_assemble(EquationSystems& es, const std::string& system_name);

    bool do_compare(ID i, ID j);

    void get_boundary_potentials(QuantumContact* qc, double& av_V, double& av_mu);

    void calculate_density(const std::string& particle);

    void init_k_space_integration(void);

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

    Device* _device;

    SimulationEnvironment* _env;

    std::set<const Boundary*> _dirichlet_boundaries;

    std::map<const Boundary*, QuantumContact*> _qc_boundaries;

    std::map<const QuantumContact*, const Boundary*> _bd_map;

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

    std::vector<unsigned int> _qc;

    NegfWrapper* _libnegf;

    options opt;

    //!diriclet DOFS
    std::set<unsigned int>  dirichlet_dofs;

    KspaceIntegration* _k_int_density;

    KspaceIntegration* _k_int_current;

    // internal status for k_space_integration
    int _which_integration;

    VectorValue<double> _k_vec;

    DofField density;

    DofField current;

    std::map<const QuantumContact*, double> _contact_potential;

    std::map<const QuantumContact*, double> _contact_current;

    double mumin, mumax;
};

#endif
