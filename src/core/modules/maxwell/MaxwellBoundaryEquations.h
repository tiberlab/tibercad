// $Id: MaxwellEquations.h 1793 2010-02-10 08:51:09Z maufder $
#ifndef _MAXWELL_BOUNDARY_EQUATIONS_H_
#define _MAXWELL_BOUNDARY_EQUATIONS_H_

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
#include "uniform_refinement_estimator.h"
#include "mesh_generation.h"
#include "mesh_modification.h"
#include "o_string_stream.h"
#include "perf_log.h"
#include "dof_map.h"
#include "numeric_vector.h"
#include "elem.h"
#include "string_to_enum.h"
#include <time.h>
#include <map>
#include "SolverException.h"
#include "SimulationInterface.h"
#include "OpticPropsInterface.h"
#include "MaxwellEquationsCommon.h"

using namespace libMesh;
using namespace std;

//! Class to solve Maxwell equations
class MaxwellBoundaryEquations : public MaxwellEquationsCommon
{
  public:
    enum Solutions
    {
      Efield,
      Efield_real,
      Efield_imag,
      Epsilon,
      Epsilon_imag,
      Mu,
      SVector,
      Intensity
    };

    //!constructor
    MaxwellBoundaryEquations(const ModelOptions& options);

    //!destructor
    virtual
    ~MaxwellBoundaryEquations(void) {};

    virtual PhysicalModel* create_physical_model(const ModelOptions& options, const Material* mat) const
            throw (ModelErrorException);

    virtual BoundaryProperties* create_boundary_model(const ModelOptions& options) const
        throw (ModelErrorException);

    static MaxwellBoundaryEquations* create(const ModelOptions& options);

    OpticPropsInterface* getOpticModel(const Elem*);

    virtual double getW() {
      return W;
    }
  protected:
    //! Get solutions at specified points in an element
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

    virtual void do_init(void);

    virtual void do_reinit(void);

    virtual void do_solve(void);

    virtual void parse_options(void) {
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

    unsigned int approxOrder;
    unsigned int extraQOrder;
    std::string inplane;// used in 2d only actually
    unsigned int defaultSourceDirection;

    std::vector<unsigned int> pmlRegions;

    double W;
    std::vector<Complex> edgeSolution;
};

inline MaxwellBoundaryEquations*
MaxwellBoundaryEquations::create(const ModelOptions& options)
{
  return (new MaxwellBoundaryEquations(options));
}

#endif
