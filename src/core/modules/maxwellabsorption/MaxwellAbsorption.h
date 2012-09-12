// $Id: MaxwellAbsorption.h 1793 2010-02-10 08:51:09Z maufder $
#ifndef _MAXWELL_ABSORPTION_H_
#define _MAXWELL_ABSORPTION_H_

#include "Device.h"
#include "mesh.h"
#include "equation_systems.h"
#include "fe.h"
#include "perf_log.h"
#include "SimulationInterface.h"
#include <map>
#include <vector>
#include <set>


/**
 * Following things are expected:
 *
 * 1. Referenced simulation is launched and solved asking solution for all integration points.
 *    Returned value is zero for each point.
 * 2. Integration points are stored.
 * 3. Now MaxwellAbsorption can be launched for solving. It launches MaxwellBoundaryEquations->solve
 *    n-times. For each integration point solution is stored.
 * 4. Referenced simulation is launched again. Now return value is real one.
 *
 */



//! Class to solve Maxwell equations
class MaxwellAbsorption : public SimulationInterface
{
  public:
    enum Solutions
    {
      Absorption
    };

    //!constructor
    MaxwellAbsorption(const ModelOptions& options);

    //!destructor
    virtual
    ~MaxwellAbsorption(void) {};

    virtual PhysicalModel* create_physical_model(const ModelOptions& options, const Material* mat) const
            throw (ModelErrorException) {
      return NULL;
    }

    virtual BoundaryProperties* create_boundary_model(const ModelOptions& options) const
        throw (ModelErrorException) {
      return NULL;
    }

    static MaxwellAbsorption* create(const ModelOptions& options);

  protected:
    SimulationInterface* maxwell;
    SimulationInterface* referencedPointsSimulation;
    std::map<Point, std::vector<Point> > integration_points; // Map which contains all integration points. elem->centroid -> element points
    std::map<Point, double> solution;

    double lambdaStart;
    double lambdaEnd;
    double lambdaStep;

    //! Get solutions at specified points in an element
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

    virtual void do_init(void);

    virtual void do_solve(void);

    virtual void parse_options(void) {
    }

  private:

    //!pointer to the device object
    static Device* _device;

    //!name of the system
    std::string system_name;

    //! Setup the available variables
    virtual void do_setup_solution_variables(void);
};

inline MaxwellAbsorption*
MaxwellAbsorption::create(const ModelOptions& options)
{
  return (new MaxwellAbsorption(options));
}

#endif
