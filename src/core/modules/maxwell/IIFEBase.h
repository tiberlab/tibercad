#ifndef __IIFEBASE_H__
#define __IIFEBASE_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "FunctionInfo.h"

class IIFEBase {
  public:
    virtual void reinit(const Elem *elem, unsigned int order, const std::vector<Point> *const pts=NULL) {
    }

    virtual void reinit(const Elem *elem, unsigned int order, const unsigned int side, const Real tolerance = TOLERANCE) {
    }

    virtual std::vector<FunctionInfo>& getFunctionsInfo() {
    }

    virtual ~IIFEBase() {
    }

    virtual void attach_quadrature_rule(QBase* q) {
    }

    virtual const std::vector<Real>& get_JxW() const {
    }
};

#endif
