#ifndef __MAXWELL_EQUATIONS_COMMON_H__
#define __MAXWELL_EQUATIONS_COMMON_H__

#include <tensor_value.h>
#include "SimulationInterface.h"

class MaxwellEquationsCommon : public SimulationInterface {
  public:
    MaxwellEquationsCommon(const ModelOptions& options)
     : SimulationInterface(options)
    {
    }

  protected:
    void addTensorSolutionR(std::vector<double>& sol, const TensorValue<Complex>& tensor, int index) {
      sol[index] = tensor(0, 0).real();
      sol[index + 1] = tensor(1, 1).real();
      sol[index + 2] = tensor(2, 2).real();
      sol[index + 3] = tensor(0, 1).real();
      sol[index + 4] = tensor(1, 2).real();
      sol[index + 5] = tensor(0, 2).real();
    }

    void addTensorSolutionI(std::vector<double>& sol, const TensorValue<Complex>& tensor, int index) {
      sol[index] = tensor(0, 0).imag();
      sol[index + 1] = tensor(1, 1).imag();
      sol[index + 2] = tensor(2, 2).imag();
      sol[index + 3] = tensor(0, 1).imag();
      sol[index + 4] = tensor(1, 2).imag();
      sol[index + 5] = tensor(0, 2).imag();
    }

  protected:

    Complex multiply(const Point& v1, const Point& v2, const TensorValue<Complex>& tensor) {

      Complex result = 0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          result += tensor(i, j) * v1(i) * v2(j);
        }
      }

      return result;
    }

    Complex multiply(const Point& v1, const VectorValue<Complex>& v2, const TensorValue<Complex> tensor) {

      Complex result = 0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          result += tensor(i, j) * v1(i) * v2(j);
        }
      }

      return result;
    }
};

#endif
