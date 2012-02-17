#ifndef __eigen_system_h__
#define __eigen_system_h__

#include "libmesh_config.h"


// Local Includes
#include "system.h"
#include "EigenSolver.h"
#include "enum_eigen_solver_type.h"
#include "EdgeDofMap.h"
#include "dense_matrix.h"
#include "IGeometryEx.h"
#include "VariableType.h"
#include <vector>

using namespace libMesh;
using namespace libMeshEnums;

/**
 * This class provides a specific system class.  It aims
 * at solving eigenvalue problems.  Currently, this class
 * is able  to handle standard eigenvalue problems
 * \p A*x=lambda*x  and generalited eigenvalue problems 
 * \p A*x=lambda*B*x.
 */

// ------------------------------------------------------------
// EigenSystem class definition

class EigenSystem : public System/*TODO , public TiberEqSystem*/
{
public:

  /**
   * Constructor.
   */
  EigenSystem (EquationSystems& es,
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
   * The type of the parent
   */
  typedef System Parent;
  
  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  virtual void clear ();

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
  virtual Complex do_get_eigen_lambda (unsigned int i);

  
  virtual Complex get_eigen_lambda (unsigned int i);

  /**
   * @returns \p "Eigen".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Eigen"; }

  virtual void get_eigen_vector(const unsigned int i, std::vector<Complex>& eigen_vector_out);

  virtual double get_relative_error(const unsigned int i);
  /**
   * @returns the number of converged eigenpairs.
   */
  unsigned int get_n_converged () const {
    return _n_converged_eigenpairs;
  }

  void setRequestedEigenPairs(unsigned int count) {
    requested_eigenpairs = count;
  }
  /**
   * @returns edge dof map
   */
  virtual EdgeDofMap* getEdgeDofMap(const bool construct = false);

  virtual void addAValue(Complex value, int i, int j) {
    matA[std::make_pair(i, j)] += value;
  }

  virtual void addBValue(Complex value, int i, int j) {
    matB[std::make_pair(i, j)] += value;
    //EigenSolver::addBRow(i, 1, &j, value);
  }

  unsigned int addVariable(const unsigned int order, const bool isVector, const unsigned int extraQOrder = 0, const bool inplane = false);

  const VariableType& getVariableType(const unsigned int i) const;

  unsigned int getVariablesCount() const;

  virtual void dof_indices(const Elem* const elem, std::vector<unsigned int>& di);

  virtual void countNonZeros(std::map<std::pair<int, int>, Complex>& values, int* result, unsigned int size);

  virtual void setSpectrumShift(const Complex& shift);

  virtual double getGeomScaling();
protected:

  virtual unsigned int initSystemSize();
  
  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * The number of converged eigenpairs.
   */
  unsigned int _n_converged_eigenpairs;

  EdgeDofMap* edgeDofMap;

  unsigned int requested_eigenpairs;

  IGeometryEx* geometryEx;

  std::map<std::pair<int, int>, Complex> matA;
  std::map<std::pair<int, int>, Complex> matB;

public:
  Complex spectrumShift;
  double solver_tolerance;
  int solver_max_it;

  //TODO Should not be here!
  SimulationInterface* simulationInterface;

  std::vector<VariableType> variables;

  IGeometryEx* getGeometryEx() {
    return geometryEx;
  }
};

#endif
