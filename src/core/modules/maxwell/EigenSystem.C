#include "EigenSystem.h"
#include "EigenSolver.h"
#include "EdgeDofMap.h"
#include "IGeometryEx.h"
#include "VariableType.h"

// ------------------------------------------------------------
// EigenSystem implementation
EigenSystem::EigenSystem (EquationSystems& es, const std::string& name, const unsigned int number) :
  System (es, name, number), _n_converged_eigenpairs (0) {
}

EigenSystem::~EigenSystem () {
}

void EigenSystem::reinit () {
  if (matA.size() > 0) {
    EigenSolver::clear_slepc();
    matA.clear();
    matB.clear();
  }
  System::reinit();
}

void EigenSystem::solve () {

  int systemSize = initSystemSize();

  std::ostringstream os;
  os << "Total system size: " << systemSize;
  Messages::debug(os.str());


  EigenSolver::prepare_slepc();

  Utils::Timer tt;

  this->assemble ();

  //get nonzeros and fill marix in EigenSolver
  int nonZerosA[systemSize];
  int nonZerosB[systemSize];

  countNonZeros(matA, nonZerosA, systemSize);
  countNonZeros(matB, nonZerosB, systemSize);

  EigenSolver::preallocate_H_matrix(systemSize, nonZerosA);

  EigenSolver::preallocate_S_matrix(systemSize, nonZerosB);

  tt.reset();

  int copyNumber = 0;
  {
    std::map<std::pair<int, int>, Complex>::iterator itA = matA.begin();
    for (; itA != matA.end(); itA++) {
      std::pair<std::pair<int, int>, Complex> pair = *itA;
      EigenSolver::addARow(pair.first.first, 1, &pair.first.second, pair.second);
      copyNumber++;
    }
    //matA.clear();
  }
  Messages::debug("Copied A"  + tt.elapsed_string());
  tt.reset();
  copyNumber = 0;

  {
    std::map<std::pair<int, int>, Complex>::iterator itB = matB.begin();
    for (; itB != matB.end(); itB++) {
      std::pair<std::pair<int, int>, Complex> pair = *itB;
      EigenSolver::addBRow(pair.first.first, 1, &pair.first.second, pair.second);
      copyNumber++;
    }
    //matB.clear();
  }
  Messages::debug("Copied B"  + tt.elapsed_string());
  //

  EigenSolver::finalize_H_assembly();
  EigenSolver::finalize_S_assembly();

  EigenSolver::SLEPCoptions options;

  options.ev_number = std::min(requested_eigenpairs, (unsigned int)(systemSize - 1));
  options.spectrum_shift = spectrumShift * getGeomScaling();
  options.eps_tolerance = solver_tolerance;
  options.matrix_output = true;

  EigenSolver::eig_value_problem_general2(options);
  
  this->_n_converged_eigenpairs = EigenSolver::number_of_converged_eigenvalues();

}

Complex EigenSystem::get_eigen_lambda (unsigned int i) {
  return do_get_eigen_lambda(i);
}

Complex EigenSystem::do_get_eigen_lambda (unsigned int i) {
  double scale = geometryEx->getScaling().get_length_scaling();

  return std::sqrt(EigenSolver::get_eigenvalue_c(i)) / scale;
}

void  EigenSystem::get_eigen_vector(const unsigned int i, std::vector<Complex>& eigen_vector_out) {
  EigenSolver::get_eigen_vector(i, eigen_vector_out);
}

void EigenSystem::setSpectrumShift(const Complex& shift) {
  spectrumShift =  shift  * shift;
}
