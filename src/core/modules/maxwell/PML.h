/*
 * PMLSEttings.h
 *
 *  Created on: May 25, 2011
 *      Author: paveryan
 */

#ifndef PML_H_
#define PML_H_

#include <vector>
#include <set>
#include "SimulationInterface.h"
#include "vector_value.h"
#include "VectorFunction.h"
#include "OpticPropsInterface.h"

class PML {
  public:
    //'Max' and 'Min' points of our device. (PML layers are not included!)
    Point minPoint;
    Point maxPoint;
    Point midPoint;

    Point allMinPoint;
    Point allMaxPoint;

    PML() {
      minPoint = Point(INT_MAX, INT_MAX);
      maxPoint = Point(INT_MIN, INT_MIN);
    }

    void init(SimulationInterface* interface) {
      MeshBase& mesh = interface->get_mesh();

      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      MeshBase::const_element_iterator el = mesh.active_local_elements_begin();


      for (; el != end_el; ++el) {
        const Elem* elem = *el;
        if (getSPML(elem, interface) <= 0) {// if it is not pml
          for (unsigned int i = 0; i < elem->n_nodes(); i++) {
            const Node node = *(elem->get_node(i));

            minPoint(0) = std::min(minPoint(0), node(0));
            minPoint(1) = std::min(minPoint(1), node(1));
            maxPoint(0) = std::max(maxPoint(0), node(0));
            maxPoint(1) = std::max(maxPoint(1), node(1));
          }
        }
      }

      midPoint = (minPoint + maxPoint) / 2;
    }

    double getSPML(const Elem* elem, SimulationInterface* interface) {
      ID subdomain = elem->subdomain_id();
      const Material* material = interface->get_environment().get_device().get_material(subdomain);

      OpticPropsInterface* opticModel =  dynamic_cast<OpticPropsInterface*>(
              material->get_model(interface->get_id()));

      return opticModel->get_spml();
    }

    bool isPMLRegion(const Elem* elem, SimulationInterface* interface) {
      return getSPML(elem, interface) > 0;
    }

    // Solution is smth like exp(ik(x+icxx/2)) ~ exp(-kcxx/2)

    libMesh::VectorValue<libMesh::Complex> getSVector(const libMesh::Point &point, double c) const {
      libMesh::VectorValue<libMesh::Complex> result(1, 1, 1);

      if (c > 0) {
        //std::cout << "ii " << point(0) << " " << minPoint(0) << " " << allMinPoint(0)<< "\n";
        for (int i = 0; i < 3; i++) {
          libMesh::Complex t, one(1, 0);
          if (point(i) > maxPoint(i)) {
            t = libMesh::Complex(1, c * (point(i) - maxPoint(i))/ (allMaxPoint(i) - maxPoint(i)) );
          } else if (point(i) < minPoint(i)) {
            t = libMesh::Complex(1, c * (point(i) - minPoint(i))/ (allMinPoint(i) - minPoint(i)) );
            // x' = x - c' * (x - xmin)^2 Where c is negative
            // In this case: c' = c/2/(xallmin - xmin), xallmin < xmin, so c must be positive
            // And we do not need '-' before c here!
          } else {
            t = 1;
          }

          result(i) =  one / t;
        }
      }

      return result;
    }

    libMesh::Complex getSVectorDet(const libMesh::Point &point, double c) const {
      libMesh::VectorValue<libMesh::Complex> sVector = getSVector(point, c);

      return sVector(0) * sVector(1) * sVector(2);
    }

    //return curls(f)
    libMesh::VectorValue<libMesh::Complex> curls(const VectorFunction &f, const libMesh::Point &point, unsigned int qp, double c) const {
      libMesh::VectorValue<libMesh::Complex> sVector = getSVector(point, c);

      return  libMesh::VectorValue<libMesh::Complex>(
          f.dphidy[qp](2)*sVector(1) - f.dphidz[qp](1)*sVector(2),
          f.dphidz[qp](0)*sVector(2) - f.dphidx[qp](2)*sVector(0),
          f.dphidx[qp](1)*sVector(0) - f.dphidy[qp](0)*sVector(1)
      );
    }
};

#endif /* PML_SETTINGS_H_ */
