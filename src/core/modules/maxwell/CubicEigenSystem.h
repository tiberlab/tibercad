#ifndef __cubic_eigen_system_h__
#define __cubic_eigen_system_h__

#include "libmesh_config.h"


// Local Includes
#include "system.h"
#include "EigenSolver.h"
#include "EigenSystem.h"
#include "enum_eigen_solver_type.h"
#include "EdgeDofMap.h"
#include "dense_matrix.h"
#include "IGeometryEx.h"
#include "VariableType.h"
#include <vector>

using namespace libMesh;
using namespace libMeshEnums;

/**
 * This class is created to solve problem: lambda*lambda*B*x = A*x + 1/(lambda1 - lambda)*C*x
 * Frist transform it to: A3x*lamba^3 + A2x*lamba^2 + A1x*lamba + A0x = 0
 */

class CubicEigenSystem : public EigenSystem {

public:

  /**
   * Constructor.
   */
    CubicEigenSystem (EquationSystems& es,
	       const std::string& name,
	       const unsigned int number);
  
    virtual Complex do_get_eigen_lambda (unsigned int i);
  
    virtual void get_eigen_vector(const unsigned int i, std::vector<Complex>& eigen_vector_out);

    // Warning global indices!
    virtual void addAValue(Complex value, int i, int j);

    virtual void addBValue(Complex value, int i, int j);

    virtual void addCValue(Complex value, int i, int j);

    virtual void setLambda1(const Complex value) {
      lambda1 = value;
    }

    virtual Complex getLambda1() {
      return lambda1;
    }
    /**
     * Assembles the system matrix.
     */
    virtual void assemble ();

    virtual void setSpectrumShift(const Complex& shift);

    virtual double getGeomScaling();
protected:

    virtual unsigned int initSystemSize();

  private:
    Complex lambda1;
    unsigned int original_system_size; // Actual is three times bigger
};

#endif
