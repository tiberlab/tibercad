#include "libmesh_config.h"

// C++ includes

// Local includes
#include "EigenSystem.h"
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
// ------------------------------------------------------------
// EigenSystem implementation
EigenSystem::EigenSystem (EquationSystems& es, const std::string& name, const unsigned int number) :
  Parent (es, name, number), _n_converged_eigenpairs (0), edgeDofMap (NULL) {
}

EigenSystem::~EigenSystem () {
  // clear data
  this->clear();
}

void EigenSystem::clear () {
  Parent::clear();
  delete edgeDofMap;
  delete geometryEx;
}

void EigenSystem::init_data () {
  Parent::init_data();
}

void EigenSystem::reinit () {
  if (matA.size() > 0) {
    EigenSolver::clear_slepc();
    matA.clear();
    matB.clear();
  }
  Parent::reinit();
}

unsigned int EigenSystem::initSystemSize () {
  return getEdgeDofMap(true)->getFreedomDegree();
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

  //TODO other options?)
  options.ev_number = requested_eigenpairs;
  options.spectrum_shift = spectrumShift.real() * getGeomScaling();
  options.eps_tolerance = solver_tolerance;

  EigenSolver::eig_value_problem_general2(options);
  
  this->_n_converged_eigenpairs = EigenSolver::number_of_converged_eigenvalues();

  std::cout << "Solved \n";
  flush(std::cout);
}

EdgeDofMap* EigenSystem::getEdgeDofMap(const bool construct) {
  if (construct || edgeDofMap == NULL) {
    if (edgeDofMap != NULL) {
      delete(edgeDofMap);
      delete geometryEx;
    }
    //TODO use init-reinit
    geometryEx = new IGeometryEx(simulationInterface, simulationInterface->get_scaling());
    //TODO bad:
    for (int i = 0; i < variables.size(); i++) {
      variables[i].scaling = simulationInterface->get_scaling().get_length_scaling();
      variables[i].dimension = get_mesh().mesh_dimension();
    }

    edgeDofMap = new EdgeDofMap(geometryEx, variables);
  }
  return edgeDofMap;
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

double EigenSystem::get_relative_error(const unsigned int i) {
  return 0.0;
}

unsigned int EigenSystem::addVariable(const unsigned int order, const bool isVector, const unsigned int extraQOrder, const bool inplane) {
  variables.push_back(VariableType(order, isVector, get_mesh().mesh_dimension(), 1.0, extraQOrder, inplane));
  return variables.size() - 1;
}

const VariableType& EigenSystem::getVariableType(const unsigned int i) const {
  return variables[i];
}

unsigned int EigenSystem::getVariablesCount() const {
  return variables.size();
}

void EigenSystem::dof_indices(const Elem* const elem, std::vector<unsigned int>& di) {
  getEdgeDofMap(false)->dof_indices(elem, di);
}

void EigenSystem::countNonZeros(std::map<std::pair<int, int>, Complex>& values, int* result, unsigned int size) {
  for (int i = 0; i < size; i++) {
    result[i] = 0;
  }

  std::map<std::pair<int, int>, Complex>::iterator it = values.begin();
  for (; it != values.end(); it++) {
    std::pair<std::pair<int, int>, Complex> pair = *it;
    result[pair.first.first]++;
  }
}

void EigenSystem::setSpectrumShift(const Complex& shift) {
  spectrumShift =  shift  * shift;
}

double EigenSystem::getGeomScaling() {
  return geometryEx->getScaling().get_length_scaling() * geometryEx->getScaling().get_length_scaling();
}

