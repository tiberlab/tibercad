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
#include "OpticPropsInterface.h"
#include "MaxwellEquationsCommon.h"
#include "ICubic.h"


using namespace std;

//! Class to solve Maxwell equations
class MaxwellEquations : public MaxwellEquationsCommon
{
  public:
    enum Solutions
    {
      EigenValue,
      EigenValue_eV,
      EigenValueImag,
      WPolariton,
      WPolariton_eV,
      WPolaritonImag,
      XHopfield,
      Epsilon,
      Epsilon_imag,
      Mu,
      SVector,
      Efield,
      Efield_real,
      Efield_imag,
      Bfield,
      Bfield_real,
      Bfield_imag,
      Poynting
      // NB: To add solution variable - add it to the start of this list
    };

    //!constructor
    MaxwellEquations(const ModelOptions& options);

    //!destructor
    virtual
    ~MaxwellEquations(void) {};

    virtual BoundaryProperties* create_boundary_model(const ModelOptions& options) const
        throw (ModelErrorException);

    static MaxwellEquations* create(const ModelOptions& options);

    virtual void plot_globaldata();

    OpticPropsInterface* getOpticModel(const Elem*);

    virtual int get_solution_for_each_mode_size() const; // Return how much solutions we have for each mode. Now it is 6: E, Ereal, Eimag, H, Himag, Hreal, Poyinting.

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
    }

    virtual void declare_E_solution(const char* name, int baseIndex, int number, bool declareOnly = false);

    virtual void declare_E_solutions(int localIndex, bool declareOnly = false);

  private:

    //!pointer to the device object
    static Device* _device;

    //!name of the system
    std::string system_name;

    //!solver options
    //options opt;

    Complex WPolaritonLow; // low polariton energy & lifetime

    std::map<ID, double> hopfieldCoeeficients; // Map: region ID -> hopfield coeff

    static void assemble_maxwell_equations(libMesh::EquationSystems& es, const std::string& system_name);

    //! Setup the available variables
    virtual void do_setup_solution_variables(void);

    //Number of requested eigen values
    unsigned int eigenCount;
    unsigned int approxOrder;
    unsigned int extraQOrder;

    Complex spectrumShift;
    double solver_tolerance;
    int solver_max_it;

    double pmlFactor;
    std::vector<bool> pmlXYZ;
    unsigned int eigensOut;
    std::string inplane;// used in 2d only actually

    bool storeSolutions;
    std::map<unsigned int, std::vector<Complex> > storedSolutions;
    
    std::vector<unsigned int> pmlRegions;
    bool polaritons;
    double lwell;
    std::string wellMaterial;
    bool useCubic;
    bool errorEstamite;

    std::map<unsigned int, unsigned int> eigenIndices;
    unsigned int accepted_eigen_count;
    bool relativeIndexing;

    void filterEigenValues();
    virtual void calculateHopfieldCoefficients();

  public:
    Complex Wc0;
    Complex Wexc0;
    double Wlt;
    std::vector<ICubic> cubics;
};

inline MaxwellEquations*
MaxwellEquations::create(const ModelOptions& options)
{
  return (new MaxwellEquations(options));
}

#endif
