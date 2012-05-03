#include "EigenSystem.h"

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

  std::cout << "Solve called \n";
  flush(std::cout);

  int systemSize = initSystemSize();

  std::cout << "Created edge dof map\n";
  flush(std::cout);
  std::cout << "Total system size: " << systemSize << "\n";
  flush(std::cout);

  EigenSolver::prepare_slepc();
  std::cout << "Prepared SLEPC \n";
  flush(std::cout);

  Utils::Timer tt;

  this->assemble ();

  std::cout << "Assemble finished " << tt.elapsed_string() << "\n";
  flush(std::cout);

  //get nonzeros and fill marix in EigenSolver
  int nonZerosA[systemSize];
  int nonZerosB[systemSize];

  countNonZeros(matA, nonZerosA, systemSize);
  countNonZeros(matB, nonZerosB, systemSize);

  std::cout << "Got nonzeros\n";
  flush(std::cout);

  EigenSolver::preallocate_H_matrix(systemSize, nonZerosA);
  std::cout << "Preallocated memory 1 \n";
  flush(std::cout);

  EigenSolver::preallocate_S_matrix(systemSize, nonZerosB);
  std::cout << "Preallocated memory 2 \n";
  flush(std::cout);

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
  std::cout << "Copied A"  << tt.elapsed_string() << " " << copyNumber << "\n";
  tt.reset();
  flush(std::cout);
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
  std::cout << "Copied B"  << tt.elapsed_string() << " " << copyNumber << "\n";
  flush(std::cout);
  //

  EigenSolver::finalize_H_assembly();
  EigenSolver::finalize_S_assembly();

  EigenSolver::SLEPCoptions options;

  options.ev_number = std::min(requested_eigenpairs, (unsigned int)(systemSize - 1));
  options.spectrum_shift = spectrumShift.real() * getGeomScaling();
  options.eps_tolerance = solver_tolerance;
  options.matrix_output = true;

  EigenSolver::eig_value_problem_general2(options);
  
  this->_n_converged_eigenpairs = EigenSolver::number_of_converged_eigenvalues();

  std::cout << "Solved \n";
  flush(std::cout);
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
