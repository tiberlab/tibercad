#include "libmesh_config.h"

// C++ includes

// Local includes
#include "EigenSystem.h"
#include "CubicEigenSystem.h"
#include "equation_systems.h"
#include "sparse_matrix.h"
#include "dof_map.h"
#include "mesh_base.h"
#include "dense_matrix.h"

#include "kelly_error_estimator.h"
#include "error_vector.h"
#include "numeric_vector.h"
#include "petsc_vector.h"
#include "EdgeDofMap.h"
#include "ElementUtils.h"
#include "MaxwellBoundaryProperties.h"
//TODO
#include "SimulationEnvironment.h"
#include "Utils.h"
// ------------------------------------------------------------
// CubicEigenSystem implementation
CubicEigenSystem::CubicEigenSystem (EquationSystems& es, const std::string& name, const unsigned int number) :
   EigenSystem (es, name, number) {
}

unsigned int CubicEigenSystem::initSystemSize () {
  original_system_size = getEdgeDofMap(true)->getFreedomDegree();
  return 3 * original_system_size;
}

void CubicEigenSystem::assemble () {
  Parent::assemble ();

  for (int i = 0; i < original_system_size; i++) {
    int i1 = i + original_system_size;
    int i2 = i + 2 * original_system_size;

    matA[std::make_pair(i1, i1)] = 1.0;
    matA[std::make_pair(i2, i2)] = 1.0;

    matB[std::make_pair(i1, i)] = 1.0;
    matB[std::make_pair(i2, i1)] = 1.0;
  }
}

Complex CubicEigenSystem::do_get_eigen_lambda (unsigned int i) {
  double scale = geometryEx->getScaling().get_length_scaling();

  return EigenSolver::get_eigenvalue_c(i) / scale;
}

void  CubicEigenSystem::get_eigen_vector(const unsigned int i, std::vector<Complex>& eigen_vector_out) {
  EigenSolver::get_eigen_vector(i, eigen_vector_out);
  eigen_vector_out.resize(original_system_size); // We do not need y & z here.
}

/**
 * This class is created to solve problem: lambda*lambda*B*x + lambda*lambda/(lambda1 - lambda)*C*x = A*x
 * Frist transform it to: A3x*lamba^3 + A2x*lamba^2 + A1x*lamba + A0x = 0.
 * A3 = B
 * A2 = -lambda1*B - C
 * A1 = -A
 * A0 = lambda1*A
 *
 * Then we solve:
 *
 * -A0  0  0   x             A1  A2  A3   x
 *  0   I  0   y  = lambda * I   0   0    y
 *  0   0  I   z             0   I   0    z
 */

void CubicEigenSystem::addAValue(Complex value, int i, int j) {
  //flush(std::cout);
  // A value goes to A0 and A1
  matA[std::make_pair(i, j)] += -lambda1 * value;
  matB[std::make_pair(i, j)] += -value;
}

void CubicEigenSystem::addBValue(Complex value, int i, int j) {
  // B value goes to A2 and A3
  matB[std::make_pair(i, j + original_system_size)] += -lambda1*value;
  matB[std::make_pair(i, j + 2*original_system_size)] += value;
}

void CubicEigenSystem::addCValue(Complex value, int i, int j) {
  // C value goes only to A2
  matB[std::make_pair(i, j + original_system_size)] += -value;
}

void CubicEigenSystem::setSpectrumShift(const Complex& shift) {
  spectrumShift =  shift;
}

double CubicEigenSystem::getGeomScaling() {
  return geometryEx->getScaling().get_length_scaling();
}
