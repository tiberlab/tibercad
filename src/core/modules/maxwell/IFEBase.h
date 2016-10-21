#ifndef __IFEBASE_H__
#define __IFEBASE_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "FunctionInfo.h"
#include "IIFEBase.h"
#include "ScalarFunction.h"
#include "ElementUtils.h"

template <class T> class IFEBase : public IIFEBase {
  protected:
  libMesh::AutoPtr<libMesh::FEBase> scalarFe;

  public://TODO
    std::vector<FunctionInfo> functions_info;
    libMesh::QBase* qBase;
    double length_scaling;
    std::vector<T> functions;
    std::vector<ScalarFunction> phi_functions;

    virtual libMesh::Point changePoint(double matrix[2][2], const libMesh::Point& point, unsigned int deriv_order) {
      libMesh::Point result = point;
      for (int i = 0; i < deriv_order; i++) {
        result = libMesh::Point(matrix[0][0] * result(0) + matrix[0][1] * result(1),
                       matrix[1][0] * result(0) + matrix[1][1] * result(1));
      }
      return result;
    }

    virtual void addFunction(const T& function, FunctionInfo::ItemType itemType, int itemId, int order) {
      functions_info.push_back(FunctionInfo(itemType, order, itemId));
      functions.push_back(function);
    }

    virtual void change_coordinate_system() {

    }

    virtual void do_reinit(const libMesh::Elem *elem, unsigned int order, const std::vector<libMesh::Point> *const pts = NULL) {
      functions_info.resize(0);
      functions.resize(0);
      phi_functions.resize(0);

      for (int i = 0; i < ElementUtils::getVertecesCount(elem); i++) {
        if (qBase != NULL || pts != NULL) {
          phi_functions.push_back(ScalarFunction((scalarFe->get_phi())[i], (scalarFe->get_dphidx())[i], (scalarFe->get_dphidy())[i], (scalarFe->get_dphidz())[i], 1.0));
        } else {
          phi_functions.push_back(ScalarFunction());
        }
      }

      std::vector<libMesh::Point> qpoints;
      if (pts != NULL) {
        qpoints = *pts;
      } else if (qBase != NULL) {
        qpoints = qBase->get_points();
      }

      addFunctions(elem, qpoints, order);
      applyScaling();
    }

    virtual void addFunctions(const libMesh::Elem *elem, const std::vector<libMesh::Point>& pts, unsigned int order) {
    }

  public:
    static const int THREE = 3; // ^^
    static const int FOUR = 4; // ^^

    IFEBase(const unsigned int dim, double scaling) : length_scaling(scaling),
        scalarFe(libMesh::FEBase::build(dim, libMesh::FEType(libMeshEnums::FIRST, libMeshEnums::LAGRANGE))), qBase(NULL) {
    }

    virtual ~IFEBase() {
      libMesh::FEBase* tmp = scalarFe.release();
      delete tmp;
    }

    virtual void attach_quadrature_rule(libMesh::QBase* q) {
      qBase = q;
      scalarFe->attach_quadrature_rule(q);
    }

    virtual const std::vector<libMesh::Real>& get_JxW() const {
      return scalarFe->get_JxW();
    }

    virtual const std::vector<libMesh::Point>& get_xyz() const {
      return scalarFe->get_xyz();
    }

    virtual void reinit(const libMesh::Elem *elem, unsigned int order, const std::vector<libMesh::Point> *const pts=NULL) {
      if (qBase != NULL || pts != NULL) {
        scalarFe->reinit(elem, pts);
      }

      do_reinit(elem, order, pts);
    }

    virtual void reinit(const libMesh::Elem *elem, unsigned int order, const unsigned int side, const libMesh::Real tolerance = libMesh::TOLERANCE) {
      scalarFe->reinit(elem, side, tolerance);

      do_reinit(elem, order);
    }

    virtual std::vector<FunctionInfo>& getFunctionsInfo() {
      return functions_info;
    }

    virtual std::vector<T>& getFunctions() {
      return functions;
    }

    virtual void applyScaling() {
    }

    virtual double getLengthScaling() const {
      return length_scaling;
    }
};

#endif
