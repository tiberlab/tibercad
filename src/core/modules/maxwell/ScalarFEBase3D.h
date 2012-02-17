#ifndef __SCALAR_FEBASE_3D_H__
#define __SCALAR_FEBASE_3D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IScalarFEBase.h"

class ScalarFEBase3D : public IScalarFEBase {
  public:
    ScalarFEBase3D(double length_scaling);

  protected:
    virtual void addTetraFunctions(const Elem *elem, const std::vector<Point>& pts, int order);

    virtual void applyScaling();

    static const unsigned int DIM = 3;

    virtual void addFunctions(const Elem *elem, const std::vector<Point>& pts, unsigned int order) {
      addTetraFunctions(elem, pts, order);
    }
};

#endif
