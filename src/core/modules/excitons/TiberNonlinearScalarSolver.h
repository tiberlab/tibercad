#ifndef _TIBERNONLINEARSCALARSOLVER_H_
#define _TIBERNONLINEARSCALARSOLVER_H_

#include <stdio.h>
#include <iostream>
#include <math.h>

/**
 * Solves f(x) = 0, where x is just a number.
 * Additionally it is supposed that f is around 1.
 * Root is in (0, 1).
 */
class TiberNonlinearScalarSolver {

  public:
    typedef void (*ScalarAssemblyRoutine)(const double& x,
                                    double* rhs,
                                    double* jacobian);

    void attach_assembly_routine(ScalarAssemblyRoutine assembly) {
      assemblyRoutine = assembly;
    }

    void rememberSolution() {
      storedSolution = solution;
    }

    void restoreSolution() {
      solution = storedSolution;
    }

    void init(double initial_solution = 0) {
      solution = initial_solution;
      rhs = 0;
      jacobian = 0;
    }

    double getSolution() {
    	return solution;
    }

    // Copy paste from TiberNonlinLS.C
    bool solve(void) {
      int max_it = 100;
      int max_step = 20;
      double rtol = 1e-7;
      double atol = 1e-20;
      double dsolution;
      double alpha;

      std::cout << "My nonlinear starts. Solution = " << solution << "\n";
      for (int k = 0; k <= max_it; k++) {
        if (k == max_it) {
          return false;
        }

        assemblyRoutine(solution, NULL, &jacobian);
        assemblyRoutine(solution, &rhs, NULL);
	std::cout << "itexc jacobian=" << jacobian << " rhs=" << rhs << "\n";

        if (fabs(rhs) < rtol * solution || fabs(rhs) < atol) {
          std::cout << "iterations: " << k << ", |du| = " << alpha * dsolution << ", |r| = " << rhs << " solution = " << solution << std::endl;
          return true;
        }

        double dsolution = rhs / jacobian;
        // the relaxation factor
        double alpha = 1.0;

        for (int ls_step = 0; ls_step <= max_step; ls_step++)
        {
          if (ls_step == max_step) {
            return false;
          }

          double new_sol = solution - alpha * dsolution;
          double new_rhs = 0;

          assemblyRoutine(new_sol, &new_rhs, NULL);

	  std::cout << "itstepexc. new_rhs= " << new_rhs << " new_solution= " << new_sol << "\n";

          if (fabs(new_rhs) < fabs(rhs)) {
            solution = new_sol;
            rhs = new_rhs;
            break;
          }

          alpha *= 0.5;
        }

        {
          std::cout << "it " << k << ", |du| = " << alpha * dsolution << ", |r| = " << rhs << "  alpha = " << alpha << " solution= " << solution << std::endl;
        }
      }

      return false;
    }

  protected:
    double storedSolution;
    double solution;
    double rhs;
    double jacobian;

    ScalarAssemblyRoutine assemblyRoutine;
};

#endif // _TIBERNONLINEARSCALARSOLVER_H_
