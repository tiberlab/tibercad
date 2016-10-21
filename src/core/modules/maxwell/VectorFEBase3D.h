#ifndef __VECTOR_FEBASE_3D_H__
#define __VECTOR_FEBASE_3D_H__

#include "fe.h"
#include "fe_type.h"
#include "quadrature.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "IVectorFEBase.h"

class VectorFEBase3D : public IVectorFEBase {
  public:
    VectorFEBase3D(double length_scaling);

  protected:
    virtual void addTetraFunctions(const libMesh::Elem *elem, const std::vector<libMesh::Point>& pts, int order);

    static const unsigned int DIM = 3;

    virtual void addFunctions(const libMesh::Elem *elem, const std::vector<libMesh::Point>& pts, unsigned int order) {
      addTetraFunctions(elem, pts, order);
    }
};

#endif
