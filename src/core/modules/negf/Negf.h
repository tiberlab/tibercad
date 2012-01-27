// $Id: Negf.h 2964 2011-10-10 20:34:57Z maufder $

#ifndef _NEGF_H_
#define _NEGF_H_

#include "SimulationInterface.h"
#include "ModelOptions.h"
#include "QuantumContact.h"
#include "TiberLinearSystem.h"
#include "Boundary.h"
#include "libnegf/NegfWrapper.h"

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
       ReorderPotential  //Laplace equation solution from reordering routines
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

    void slice(void);

  protected:

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    //virtual void do_setup_solution_variables(void);


    //! Solve the MyPoisson equation
    virtual void do_solve(void);

    void do_assemble(const ModelOptions& opt);

    void print_ham(std::string form);

    void print_Lib(void);

    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
                                                  const Material* mat) const;

    //! We need to create boundary condition model
    //virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
    //    const MaterialBoundary* boundary) const;

    virtual void do_setup_solution_variables(void);


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

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

    void get_blocks(std::vector<int>& cblok);

    struct options
    {
        unsigned int n_blocks;

        std::string pot_module;

        double Emin;

        double Emax;

        double Estep;

        std::vector <int> Np_n;

        double Np_real;


    };



  private:

    void activate_quantum_contacts(void);
    void deactivate_quantum_contacts(void);

    Device* _device;

    std::map<const Boundary*, ID> _boundaries;

    SimulationEnvironment* _env;

    std::map<ID, QuantumContact*> _quantum_contacts;

    Negf(const ModelOptions& option);

    TiberLinearSystem* _sys;
    TiberLinearSystem* _sys_H;
    TiberLinearSystem* _sys_S;

    std::vector<ID> _permu;

    std::vector<ID> _end_blocks;

    unsigned int _device_n_dofs;

    std::vector<unsigned int> _dev;

    std::vector<unsigned int> _qc_n_dofs;

    std::vector<unsigned int> _qc;

    NegfWrapper* _libnegf;

    options opt;

};

#endif
