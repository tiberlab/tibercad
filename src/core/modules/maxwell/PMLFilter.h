/*
 * PMLSEttings.h
 *
 *  Created on: May 25, 2011
 *      Author: paveryan
 */

#ifndef PML_FILTER_H_
#define PML_FILTER_H_

#include <vector>
#include <set>
#include "SimulationInterface.h"
#include "MaxwellEquations.h"
#include "EigenSolver.h"

class PMLFilter {
  private:
    double factor;
    SimulationInterface* maxwell;

  public:
    PMLFilter(double value, SimulationInterface* maxwellEquations) {
      factor = value;
      maxwell = maxwellEquations;
    }

    void filter(EigenSystem& system, std::map<unsigned int, unsigned int>& eigenIndices, std::map<unsigned int, std::vector<Complex> >& storedSolutions, std::vector<bool>& pmlXYZ) {
      eigenIndices.clear();

      std::map<ID, std::vector<double> > solutions;

      for (int i = 0; i < system.get_n_converged(); i++) {
        solutions[MaxwellEquations::Efield + 3*i].resize(0);
      }

      MeshBase& mesh = maxwell->get_mesh();

      MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

      PML& pml = system.getGeometryEx()->pml;

      std::vector<double> maxValue(system.get_n_converged(), 0);
      std::vector<double> maxBorderPMLValue(system.get_n_converged(), 0);


      for (; el != end_el; ++el) {
        const Elem* elem = *el;

        maxwell->get_solution(elem, solutions);

        bool isPML = pml.isPMLRegion(elem, maxwell); // respect to +, * etc. slow operation
        for (int i = 0; i < system.get_n_converged(); i++) {
          for (int qp = 0; qp < elem->n_nodes(); qp++) {
            std::vector<double>& solution = solutions[MaxwellEquations::Efield + 3*i];
            double E = std::sqrt(solution[3*qp] * solution[3*qp] +
                                 solution[3*qp + 1] * solution[3*qp + 1] +
                                 solution[3*qp + 2] * solution[3*qp + 2]);
            if (isPML && isOnPmlBorder(*(elem->get_node(qp)), pml, pmlXYZ)) {
              maxBorderPMLValue[i] = std::max(maxBorderPMLValue[i], E);
            } else {
              maxValue[i] = std::max(maxValue[i], E);
            }
          }
        }
      }

      for (int i = system.get_n_converged() - 1; i >= 0; i--) {
        bool acceptValue = true;

        Complex eigen = system.get_eigen_lambda(i);

        acceptValue = eigen.real() > 0;

        if (acceptValue) {
          acceptValue = maxBorderPMLValue[i] < maxValue[i] * factor;
        }

        if (acceptValue) {
          unsigned int num = eigenIndices.size();
          eigenIndices.insert(std::make_pair(num, i));
        } else {
          storedSolutions.erase(i); //Discard this solution.
        }
      }
    }

    bool isOnPmlBorder(const Node& node, PML& pml, std::vector<bool>& pmlXYZ) {
      for (int i = 0; i < 3; i++) {
        if (!pmlXYZ[i]) {
          return false;
        }

        if (std::abs((node(i) - pml.allMinPoint(i))/(pml.minPoint(i) - pml.allMinPoint(i))) < 0.5 ||
            std::abs((node(i) - pml.allMaxPoint(i))/(pml.maxPoint(i) - pml.allMaxPoint(i))) < 0.5) {
          return true;
        }
      }
      return false;
    }
};

#endif /* PML_FILTER_H_ */
