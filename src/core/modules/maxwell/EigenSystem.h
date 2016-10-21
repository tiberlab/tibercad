#ifndef __eigen_system_h__
#define __eigen_system_h__

#include "libmesh_config.h"


// Local Includes
#include "VectorSystem.h"
#include <vector>


using namespace libMeshEnums;

namespace libMesh{
class EquationSystems;
}

/**
 * This class provides a specific system class.  It aims
 * at solving eigenvalue problems.  Currently, this class
 * is able  to handle standard eigenvalue problems
 * \p A*x=lambda*x  and generalited eigenvalue problems 
 * \p A*x=lambda*B*x.
 */

// ------------------------------------------------------------
// EigenSystem class definition

class EigenSystem : public VectorSystem, public libMesh::System {

  public:

    /**
     * Constructor.
     */
    EigenSystem (libMesh::EquationSystems& es,
	       const std::string& name,
	       const unsigned int number);

    /**
     * Destructor.
     */
    virtual ~EigenSystem ();
 
    /**
     * The type of system.
     */
    typedef EigenSystem sys_type;

    /**
     * @returns a clever pointer to the system.
     */
    sys_type & system () { return *this; }

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
   * Assembles & solves the eigen system.
   */

  /**
   * Returns real and imaginary part of the ith eigenvalue and copies
   * the respective eigen vector to the solution vector.
   */
  virtual libMesh::Complex do_get_eigen_lambda (unsigned int i);

  
  virtual libMesh::Complex get_eigen_lambda (unsigned int i);

  /**
   * @returns \p "Eigen".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "VectorEigen"; }

  virtual void get_eigen_vector(const unsigned int i, std::vector<libMesh::Complex>& eigen_vector_out);

  /**
   * @returns the number of converged eigenpairs.
   */
  unsigned int get_n_converged () const {
    return _n_converged_eigenpairs;
  }

  void setRequestedEigenPairs(unsigned int count) {
    requested_eigenpairs = count;
  }

  virtual void addAValue(libMesh::Complex value, int i, int j) {
    matA[std::make_pair(i, j)] += value;
  }

  virtual void addBValue(libMesh::Complex value, int i, int j) {
    matB[std::make_pair(i, j)] += value;
    //EigenSolver::addBRow(i, 1, &j, value);
  }

  virtual void setSpectrumShift(const libMesh::Complex& shift);

protected:

  /**
   * The number of converged eigenpairs.
   */
  unsigned int _n_converged_eigenpairs;

  unsigned int requested_eigenpairs;

  std::map<std::pair<int, int>, libMesh::Complex> matA;
  std::map<std::pair<int, int>, libMesh::Complex> matB;

public:
  libMesh::Complex spectrumShift;
  double solver_tolerance;
  int solver_max_it;
};

#endif
