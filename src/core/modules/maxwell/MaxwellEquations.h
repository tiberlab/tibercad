// $Id: MaxwellEquations.h 1793 2010-02-10 08:51:09Z maufder $
#ifndef _MAXWELLEQUATIONS_H_
#define _MAXWELLEQUATIONS_H_

#include "Device.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "sparse_matrix.h"
#include "mesh_refinement.h"
#include "error_vector.h"
#include "patch_recovery_error_estimator.h"
#include "uniform_refinement_estimator.h"
#include "hp_coarsentest.h"
#include "hp_singular.h"
#include "mesh_generation.h"
#include "mesh_modification.h"
#include "o_string_stream.h"
#include "perf_log.h"
#include "getpot.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "elem.h"
#include "string_to_enum.h"
#include "gmsh_io.h"
#include "EigenSolver.h"
#include "CubicEigenSystem.h"
#include <time.h>
#include <map>
#include "ElementUtils.h"
#include "SolverException.h"
#include "OpticPropsModel.h"


using namespace libMesh;
using namespace std;

//! Class to solve Maxwell equations
class MaxwellEquations : public SimulationInterface
{
  public:
    enum Solutions
    {
      EigenValue,
      EigenValueImag,
      XHopfield,
      Epsilon,
      Epsilon_imag,
      Mu,
      SVector,
      Efield // MUST BE LAST
    };

    //!constructor
    MaxwellEquations(const ModelOptions& options);

    //!destructor
    virtual
    ~MaxwellEquations(void) {};

    virtual PhysicalModel* create_physical_model(const ModelOptions& options, const Material* mat) const
            throw (ModelErrorException);

    virtual BoundaryProperties* create_boundary_model(const ModelOptions& options) const
        throw (ModelErrorException);

    static MaxwellEquations* create(const ModelOptions& options);

    virtual void plot_globaldata();

    OpticPropsModel* getOpticModel(const Elem*);
  protected:
    //virtual void build_nodal_results(const std::set<std::string>& variables, std::vector<
        //double>& results, std::vector<std::string>& legend);
    //! Get solutions not located on the mesh or the atoms
    virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);

    //! Get solutions at specified points in an element
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

    virtual void do_init(void);

    virtual void do_solve(void);

    virtual void parse_options(void) {
      std::cout << "EEEPO"; flush(std::cout);
    }

  private:

    //!pointer to the device object
    static Device* _device;

    //!name of the system
    std::string system_name;

    //!solver options
    //options opt;

    static void assemble_maxwell_equations(EquationSystems& es, const std::string& system_name);

    //! Setup the available variables
    virtual void do_setup_solution_variables(void);

    //Number of requested eigen values
    unsigned int eigenCount;
    unsigned int approxOrder;
    unsigned int extraQOrder;

    double spectrumShift;
    double solver_tolerance;
    int solver_max_it;

    unsigned int maxIterations;
    double relativeError;
    unsigned int eigensOut;
    std::string inplane;// used in 2d only actually

    bool storeSolutions;
    std::map<unsigned int, std::vector<Complex> > storedSolutions;
    
    std::map<Point, double> hopfieldCoefficients;
    std::vector<unsigned int> pmlRegions;
    bool polaritons;

    std::map<unsigned int, unsigned int> eigenIndices;
    unsigned int accepted_eigen_count;
    bool relativeIndexing;

    void filterEigenValues(double factor);

  public:
    Complex Wc;
    Complex Wexc0;
    double Wlt;
    Complex VRabiApprox;

    virtual double getXHopfield(int i);
};

inline MaxwellEquations*
MaxwellEquations::create(const ModelOptions& options)
{
  return (new MaxwellEquations(options));
}

#endif
