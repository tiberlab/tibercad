#ifndef __VARIABLE_TYPE_H_
#define __VARIABLE_TYPE_H_

#include "IIFEBase.h"
#include "VectorFEBase1D.h"
#include "VectorFEBase2D.h"
#include "VectorFEBase3D.h"
#include "ScalarFEBase1D.h"
#include "ScalarFEBase2D.h"
#include "ScalarFEBase3D.h"
#include "VectorFEBase2DOutplane.h"

class VariableType {
  public:
    unsigned int order;

    //IIFEBase* feBase;

    unsigned int dimension;

    bool isVector;

    double scaling;

    unsigned int extraQOrder;

    bool inplane;

    VariableType(unsigned int o, bool isV, unsigned int dim, double length_scaling, unsigned int extraQO, bool flag): order(o), isVector(isV), scaling(length_scaling), dimension(dim), extraQOrder(extraQO), inplane(flag) {
    }

    ~VariableType() {
      //delete feBase; TODO
    }

    IIFEBase* getFEbase() const {
      if (isVector && dimension == 2) {
        if (inplane) {
          return new VectorFEBase2D(scaling);
        } else {
          return new VectorFEBase2DOutplane(scaling);
        }
      } else if (!isVector && dimension == 2) {
        return new ScalarFEBase2D(scaling);
      } else if (dimension == 1) {
        return new VectorFEBase1D(scaling);
      } else if (isVector) {
        return new VectorFEBase3D(scaling);
      } else {
        return new ScalarFEBase3D(scaling);
      }
    }
};

#endif /* __VARIABLE_TYPE_H_ */

