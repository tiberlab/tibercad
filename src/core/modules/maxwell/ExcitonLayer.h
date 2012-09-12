/*
 * ExcitonLayer.h
 *
 *  Created on: Sep 26, 2011
 *      Author: paveryan
 */

#ifndef EXCITONLAYER_H_
#define EXCITONLAYER_H_

#include "EigenSystem.h"
#include "IGeometryEx.h"
#include "MaxwellEquations.h"
#include "Constants.h"
#include "PML.h"
#include "OpticPropsInterface.h"

class ExcitonLayer
{
  public:
    Complex Wexc;
    double Q;// This is Q/c in fact
    MaxwellEquations* mme;

    ExcitonLayer(double length_scaling, MaxwellEquations* mm) {
      mme = mm;

      Wexc = mme->Wexc0;
      Q = mme->Wlt; // real value of Q TODO: add well size and Bohr exciton radius

      //Scaling
      Wexc = (Wexc / Constants::c) * length_scaling;
      Q = (Q / Constants::c) * length_scaling;

    }

    virtual ~ExcitonLayer() {
    }

    void preset(EigenSystem& system) {
      CubicEigenSystem& cubicSystem = dynamic_cast<CubicEigenSystem&>(system);
      cubicSystem.setLambda1(Wexc);
    }

    void addData(EigenSystem& system) {
      // find the exciton simulation to use
      SimulationInterface* _exciton_sim = mme->find_simulation("excitontransport");

      if (_exciton_sim == NULL) {
        std::string msg("ExcitonLayer: Simulation not found");
        throw InitFailedException(msg);
      } else if (!_exciton_sim->is_solved() || mme->Wc0.real() == -1) {
        return;
      }

      ID densityId = _exciton_sim->get_solution_id("Xdens");

      CubicEigenSystem& cubicSystem = dynamic_cast<CubicEigenSystem&>(system);

      std::map<unsigned int, Point> FexcFi;
      std::map<unsigned int, Point> FexcFiEpsilon;
      double Fexc = 0.0;

      const MeshBase& mesh = system.get_mesh();
      const unsigned int dimension = mesh.mesh_dimension();

      IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
      const VariableType fe_type = system.getVariableType(0);
      //TODO to be revised
      AutoPtr<QBase> qrule(new QGauss(dimension, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

      MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

      for (; el != end_el; ++el) {
        const Elem* elem = *el;

        PML pml = system.getGeometryEx()->pml;

        OpticPropsInterface* opticModel = mme->getOpticModel(elem);

        if (!pml.isPMLRegion(elem, mme)) {

          std::vector<unsigned int> all_dof_indices;
          system.dof_indices (elem, all_dof_indices);

          fe->attach_quadrature_rule (qrule.get());
          fe->reinit (elem, system.getEdgeDofMap(false)->getPOrder(elem, 0));

          const std::vector<Real>& JxW = fe->get_JxW();
          const std::vector<VectorFunction >& edge_phi = fe->getFunctions();
          const std::vector<Point>& xyz = fe->get_xyz();

          for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
            const Point& point = xyz[qp];
            double Xden = 0.0;//std::abs(point(0)) < 2 ? 1 : 0;
            _exciton_sim->get_solution(elem, densityId, Xden, point);

            Fexc += JxW[qp] * Xden;

            for (unsigned int i=0; i<edge_phi.size(); i++) {
              if (all_dof_indices[i] != ElementUtils::INVALID_FUNCTION_ID) {
                 Point value = JxW[qp] * std::sqrt(Xden) * edge_phi[i].phi[qp];

                 if (Xden != 0) {
                   if (FexcFi.find(all_dof_indices[i]) == FexcFi.end()) {
                     FexcFi[all_dof_indices[i]] = value;
                     FexcFiEpsilon[all_dof_indices[i]] = value * opticModel->get_dielectric_constant().real();//TODO complex?
                   } else {
                     FexcFi[all_dof_indices[i]] += value;
                     FexcFiEpsilon[all_dof_indices[i]] += value * opticModel->get_dielectric_constant().real();//TODO complex?
                   }
                 }
              }
            }
          }
        }
      }

      delete fe;

      std::map<unsigned int, Point>::iterator it1;
      std::map<unsigned int, Point>::iterator it2;

      for (it1 = FexcFi.begin(); it1 != FexcFi.end(); it1++) {
        for (it2 = FexcFiEpsilon.begin(); it2 != FexcFiEpsilon.end(); it2++) {
          Point p1 = (*it1).second;
          Point p2 = (*it2).second;

          cubicSystem.addCValue(Q * (p1 * p2) / Fexc, (*it1).first, (*it2).first);
        }
      }
    }

    // Here epsilon is considered to be constant

    Complex getApproximateRabiSplitting(EigenSystem& system) {
      SimulationInterface* _exciton_sim = mme->find_simulation("excitontransport");

      if (_exciton_sim == NULL) {
        std::string msg("ExcitonLayer: Simulation not found");
        throw InitFailedException(msg);
      } else if (!_exciton_sim->is_solved() || mme->Wc0.real() == -1) {
        return -100;
      }

      ID densityId = _exciton_sim->get_solution_id("Xdens");
      ID E0id = mme->get_solution_id("Efield_0");

      double fi0I = 0.0; // integral of sqrt(Xdens)*E
      double E2I = 0.0;  // integral of E*E
      double XdensI = 0.0; // integral of Xdens



      const MeshBase& mesh = system.get_mesh();
      const unsigned int dimension = mesh.mesh_dimension();

      IVectorFEBase* fe = dynamic_cast<IVectorFEBase*>(system.getVariableType(0).getFEbase());
      const VariableType fe_type = system.getVariableType(0);

      AutoPtr<QBase> qrule(new QGauss(dimension, static_cast<Order>(2 * fe_type.order + 2 + fe_type.extraQOrder)));

      MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

      for (; el != end_el; ++el) {
        const Elem* elem = *el;

        std::vector<unsigned int> all_dof_indices;
        system.dof_indices (elem, all_dof_indices);

        fe->attach_quadrature_rule (qrule.get());
        fe->reinit (elem, system.getEdgeDofMap(false)->getPOrder(elem, 0));

        const std::vector<Real>& JxW = fe->get_JxW();
        const std::vector<VectorFunction >& edge_phi = fe->getFunctions();
        const std::vector<Point>& xyz = fe->get_xyz();

        double E[3];

        std::vector<double> Esol;

        mme->get_solution(elem, E0id, Esol, xyz);

        for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
          const Point& point = xyz[qp];
          double Xden = 0.0;//std::abs(point(0)) < 2 ? 1 : 0;
          _exciton_sim->get_solution(elem, densityId, Xden, point);

          XdensI += JxW[qp] * std::sqrt(Xden) * std::sqrt(Xden);
          E2I += JxW[qp] * Esol[3*qp + 1] * Esol[3*qp + 1]; // TODO TODO

          fi0I += JxW[qp] * Esol[3*qp + 1] * std::sqrt(Xden);
        }
      }

      delete fe;

      std::cout << "PPAPP " << std::sqrt(mme->Wlt * mme->Wexc0 * (fi0I * fi0I / E2I / XdensI / 2)) << "\n";
      return std::sqrt(mme->Wlt * mme->Wexc0 * (fi0I * fi0I / E2I / XdensI / 2));

    }
};

#endif /* EXCITONLAYER_H_ */
