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

class PMLFilter {
  private:
    double factor;
    SimulationInterface* maxwell;

  public:
    PMLFilter(double value, SimulationInterface* maxwellEquations) {
      factor = value;
      maxwell = maxwellEquations;
    }

    void filter(EigenSystem& system, std::map<unsigned int, unsigned int>& eigenIndices, std::map<unsigned int, std::vector<Complex> >& storedSolutions) {
      eigenIndices.clear();

      std::map<ID, std::vector<double> > solutions;

      for (int i = 0; i < system.get_n_converged(); i++) {
        solutions[MaxwellEquations::Efield + i].resize(0);
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
            std::vector<double>& solution = solutions[MaxwellEquations::Efield + i];
            double E = std::sqrt(solution[3*qp] * solution[3*qp] +
                                 solution[3*qp + 1] * solution[3*qp + 1] +
                                 solution[3*qp + 2] * solution[3*qp + 2]);
            if (isPML && isOnPmlBorder(*(elem->get_node(qp)), pml)) {
              maxBorderPMLValue[i] = std::max(maxBorderPMLValue[i], E);
            } else {
              maxValue[i] = std::max(maxValue[i], E);
            }
          }
        }
      }

      for (int i = system.get_n_converged() - 1; i >= 0; i--) {
        bool acceptValue = false;

        Complex eigen = system.get_eigen_lambda(i);

        acceptValue = eigen.real() > 0;

        if (acceptValue) {
          acceptValue = maxBorderPMLValue[i] < maxValue[i] * 0.1;
        }

        if (acceptValue) {
          unsigned int num = eigenIndices.size();
          eigenIndices.insert(std::make_pair(num, i));
        } else {
          storedSolutions.erase(i); //Discard this solution.
        }
      }
    }

    bool isOnPmlBorder(const Node& node, PML& pml) {
      for (int i = 0; i < 3; i++) {
        if (std::abs((node(i) - pml.allMinPoint(0))/(pml.minPoint(0) - pml.allMinPoint(0))) < 0.5 ||
            std::abs((node(i) - pml.allMaxPoint(0))/(pml.maxPoint(0) - pml.allMaxPoint(0))) < 0.5) {
          return true;
        }
      }
      return false;
    }
};

#endif /* PML_FILTER_H_ */
