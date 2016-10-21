#ifndef __VECTOR_LINEAR_SYSTEM_h__
#define __VECTOR_LINEAR_SYSTEM_h__

#include "libmesh_config.h"


// Local Includes
#include "VectorSystem.h"
#include "EdgeDofMap.h"
#include "IGeometryEx.h"
#include "VariableType.h"
#include <vector>
#include "TiberLinearSystem.h"
#include "linear_implicit_system.h"


using namespace libMeshEnums;

/**
 * This class provides a specific system class.  It aims
 * at solving linear problems Ax=b, gives possibility to use custom DofMap object for vector finite elements.
 * Also solves complex problem.
 * Lack of multiprocess libMesh staff.
 */

class VectorLinearSystem : public VectorSystem, public TiberLinearSystem {
  public:

  /**
   * Constructor.
   */
    VectorLinearSystem (libMesh::EquationSystems& es,
	       const std::string& name,
	       const unsigned int number);

    /**
     * Destructor.
     */
    virtual ~VectorLinearSystem ();

    /**
     * @returns a clever pointer to the system.
     */
    VectorLinearSystem & system () { return *this; }

    /**
     * Reinitializes the member data fields associated with
     * the system, so that, e.g., \p assemble() may be used.
     */
    virtual void reinit ();

    /**
     * Assembles & solves the eigen system.
     */
    virtual void solve ();

    /**
     * @returns \p "VectorLinear".  Helps in identifying
     * the system type in an equation system file.
     */
    virtual std::string system_type () const { return "VectorLinear"; }

    virtual void addAValue(libMesh::Complex value, int i, int j) {
      matA[std::make_pair(i, j)] += value;
    }

    virtual void addBValue(libMesh::Complex value, int i) {
      columnB[i] += value;
    }

    virtual void get_solution(std::vector<libMesh::Complex>& result);
  protected:
    virtual void init_data();

    std::map<std::pair<int, int>, libMesh::Complex> matA;
    std::map<int, libMesh::Complex> columnB;


  private:
    unsigned int systemSize;
};

#endif
