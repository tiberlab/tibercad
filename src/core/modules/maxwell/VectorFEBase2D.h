#ifndef __VECTOR_FEBASE_2D_H__
#define __VECTOR_FEBASE_2D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "VectorFunction.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IVectorFEBase.h"

class VectorFEBase2D : public IVectorFEBase {
  public:
    VectorFEBase2D(double length_scaling);

  protected:
    virtual void change_coordinate_system();

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
