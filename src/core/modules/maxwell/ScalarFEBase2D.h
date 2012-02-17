#ifndef __SCALAR_FEBASE_2D_H__
#define __SCALAR_FEBASE_2D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IScalarFEBase.h"

class ScalarFEBase2D : public IScalarFEBase {
  public:
    ScalarFEBase2D(double length_scaling);

  protected:
    void change_coordinate_system();

    virtual void addTriFunctions(const Elem *elem, const std::vector<Point>& pts, int order);

    virtual void addQuadFunctions(const Elem *elem, const std::vector<Point>& pts, int order);

    static const unsigned int THREE = 3; // ^^

    static const unsigned int FOUR = 4; // ^^

    static const unsigned int DIM = 2;

    virtual void addFunctions(const Elem *elem, const std::vector<Point>& pts, unsigned int order) {
      if (ElementUtils::getVertecesCount(elem) == THREE) {
        addTriFunctions(elem, pts, order);
      } else {
        addQuadFunctions(elem, pts, order);
      }
    }
};

#endif
