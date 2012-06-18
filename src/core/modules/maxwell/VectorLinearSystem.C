#include "VectorLinearSystem.h"
#include <sparse_matrix.h>
#include <linear_solver.h>
#include "TiberLinearSolver.h"

// ------------------------------------------------------------
// EigenSystem implementation
VectorLinearSystem::VectorLinearSystem (EquationSystems& es, const std::string& name, const unsigned int number) :
  TiberLinearSystem (es, name, number) {
}

VectorLinearSystem::~VectorLinearSystem () {
}

void VectorLinearSystem::reinit () {
  //TiberLinearSystem::reinit();
}

void VectorLinearSystem::solve () {
  systemSize = initSystemSize();

  // Resize the solution conformal to the current mesh
  solution->init (2*systemSize, 2*systemSize, false, PARALLEL);

  if (matrix == NULL) {
    matrix = &(this->add_matrix ("System Matrix"));
  }

  matrix->init(2*systemSize, 2*systemSize, 2*systemSize, 2*systemSize);

  if (rhs == NULL)
    rhs = &(this->add_vector ("RHS Vector", false));

  rhs->init (2*systemSize, 2*systemSize, false, PARALLEL);


  /***************************************************************************/
  this->assemble ();
  // TODO for one test
  //addBValue(Complex(1, 0), (systemSize - 1) / 2);

  // We solve Ax = b
  // (A1+iA2)(x1 +ix2) = b1 + ib2
  // A1x1-A2x2 = b1
  // A1x2+A2x1 = b2
  {
    std::map<std::pair<int, int>, Complex>::iterator itA = matA.begin();
    for (; itA != matA.end(); itA++) {
      std::pair<std::pair<int, int>, Complex> pair = *itA;
      matrix->set(pair.first.first, pair.first.second, pair.second.real());
      matrix->set(pair.first.first, pair.first.second + systemSize, - pair.second.imag());

      matrix->set(systemSize + pair.first.first, pair.first.second, pair.second.imag());
      matrix->set(systemSize + pair.first.first, systemSize + pair.first.second, pair.second.real());

      //std::cout << "Matrix output " << pair.first.first << " " << pair.first.second << " " << (pair.second) << "\n";
      //std::cout << "Matrix output " << (systemSize + pair.first.first) << " " << (systemSize + pair.first.second) << " " << (pair.second.real() + pair.second.imag()) << "\n";
    }
  }

  {
    std::map<int, Complex>::iterator itB = columnB.begin();
    for (; itB != columnB.end(); itB++) {
      std::pair<int, Complex> pair = *itB;
      rhs->set(pair.first, pair.second.real());
      rhs->set(systemSize + pair.first, pair.second.imag());
      //std::cout << "Matrix output B" << pair.first << " " << pair.second << " " << "\n";
    }
  }

/*
  std::pair<unsigned int, Real> result = linear_solver->solve (*matrix, this->request_matrix("Preconditioner"), *solution, *rhs, 1e-30, 1000);



  std::cout << "RETURN INFO " << result.first << " " << result.second << "\n";
*/

/*
  std::cout << "RHS1 " << rhs->l1_norm() << " " << rhs->l2_norm() << "\n";
  matrix->vector_mult(*rhs, *solution);
  std::cout << "RHS2 " << rhs->l1_norm() << " " << rhs->l2_norm() << "\n";
*/

/*
  NumericVector<Number>* rhs2 = &(this->add_vector ("RHS Vector2", false));
  rhs2->init (2*systemSize, 2*systemSize, false, PARALLEL);
  {
    std::map<int, Complex>::iterator itB = columnB.begin();
    for (; itB != columnB.end(); itB++) {
      std::pair<int, Complex> pair = *itB;
      rhs2->set(pair.first, - pair.second.real());
      rhs2->set(systemSize + pair.first, - pair.second.imag());
      //std::cout << "Matrix output B" << pair.first << " " << pair.second << " " << "\n";
    }
  }
*/

  //std::cout << typeid(*linear_solver).name() << "\n";

/*
  std::cout << "RHS3 " << rhs2->l1_norm() << " " << rhs2->l2_norm() << "\n";
  matrix->vector_mult_add(*rhs2, *solution);
  std::cout << "RHS4 " << rhs2->l1_norm() << " " << rhs2->l2_norm() << "\n";
*/
  ModelOptions options = ModelOptions();
  options.set_option("relative_tolerance", 1e-12);
  options.set_option("preconditioner", "lu");
  //options.


  TiberLinearSolver* tiber_solver = TiberLinearSolver::create(options);
  tiber_solver->init();
  std::pair<unsigned int, Real> result2 = tiber_solver->solve(*matrix, *(solution.get()), *rhs);
  tiber_solver->clear();
  delete tiber_solver;

  //std::cout << "RETURN INFO 2" << result2.first << " " << result2.second << "\n";

}

void VectorLinearSystem::init_data() {
}

void VectorLinearSystem::get_solution(std::vector<Complex>& result) {
  result.resize(systemSize);
  for (int i = 0; i < systemSize; i++) {
    Complex tmp(solution->el(i), solution->el(systemSize + i));
    result[i] = tmp;
    //std::cout << "Solution " << i << " : " << tmp << "\n";
  }
}
