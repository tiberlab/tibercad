#ifndef __IIFEBASE_H__
#define __IIFEBASE_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "FunctionInfo.h"

class IIFEBase {
  public:
    virtual void reinit(const libMesh::Elem *elem, unsigned int order, const std::vector<libMesh::Point> *const pts=NULL) {
    }

    virtual void reinit(const libMesh::Elem *elem, unsigned int order, const unsigned int side, const libMesh::Real tolerance = libMesh::TOLERANCE) {
    }

    virtual std::vector<FunctionInfo>& getFunctionsInfo() {
    }

    virtual ~IIFEBase() {
    }

    virtual void attach_quadrature_rule(libMesh::QBase* q) {
    }

    virtual const std::vector<libMesh::Real>& get_JxW() const {
    }
};

#endif
