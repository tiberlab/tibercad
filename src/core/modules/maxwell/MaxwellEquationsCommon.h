#ifndef __MAXWELL_EQUATIONS_COMMON_H__
#define __MAXWELL_EQUATIONS_COMMON_H__

#include <tensor_value.h>
#include "SimulationInterface.h"
#include "OpticPropsModel.h"
#include "OpticPropsInterface.h"

class MaxwellEquationsCommon : public SimulationInterface {
  public:
    MaxwellEquationsCommon(const ModelOptions& options)
     : SimulationInterface(options)
    {
    }

    virtual PhysicalModel* create_physical_model(const ModelOptions& options, const Material* mat) const
            throw (ModelErrorException) {
      const std::string& modelname = get_option((mat->get_name() + "_opticmodel"), "");

      if (modelname == "") {
        return OpticPropsModel::create(options);
      } else {
        OpticPropsInterface* model =
            OpticPropsInterface::create(modelname, mat, options);

        if (model == NULL)
          throw ModelErrorException(
              "Maxwell: No such physical model: " + modelname);

        return model;
      }
    }


    static void addVectorSolutionR(std::map<ID, std::vector<double> >& sol, ID id, const VectorValue<Complex>& vector, int index) {
      if (sol.count(id)) {
        sol[id][index] = vector(0).real();
        sol[id][index + 1] = vector(1).real();
        sol[id][index + 2] = vector(2).real();
      }
    }

    static void addVectorSolutionI(std::map<ID, std::vector<double> >& sol, ID id, const VectorValue<Complex>& vector, int index) {
      if (sol.count(id)) {
        sol[id][index] = vector(0).imag();
        sol[id][index + 1] = vector(1).imag();
        sol[id][index + 2] = vector(2).imag();
      }
    }

    static void addVectorSolutionA(std::map<ID, std::vector<double> >& sol, ID id, const VectorValue<Complex>& vector, int index) {
      if (sol.count(id)) {
        sol[id][index] = std::abs(vector(0));
        sol[id][index + 1] = std::abs(vector(1));
        sol[id][index + 2] = std::abs(vector(2));
      }
    }

    static void addTensorSolutionR(std::map<ID, std::vector<double> >& sol, ID id, const TensorValue<Complex>& tensor, int index) {
      if (sol.count(id)) {
        sol[id][index] = tensor(0, 0).real();
        sol[id][index + 1] = tensor(1, 1).real();
        sol[id][index + 2] = tensor(2, 2).real();
        sol[id][index + 3] = tensor(0, 1).real();
        sol[id][index + 4] = tensor(1, 2).real();
        sol[id][index + 5] = tensor(0, 2).real();
      }
    }

    static void addTensorSolutionI(std::map<ID, std::vector<double> >& sol, ID id, const TensorValue<Complex>& tensor, int index) {
      if (sol.count(id)) {
        sol[id][index] = tensor(0, 0).imag();
        sol[id][index + 1] = tensor(1, 1).imag();
        sol[id][index + 2] = tensor(2, 2).imag();
        sol[id][index + 3] = tensor(0, 1).imag();
        sol[id][index + 4] = tensor(1, 2).imag();
        sol[id][index + 5] = tensor(0, 2).imag();
      }
    }

    static Complex multiply(const Point& v1, const Point& v2, const TensorValue<Complex>& tensor) {

      Complex result = 0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          result += tensor(i, j) * v1(i) * v2(j);
        }
      }

      return result;
    }

    static Complex multiply(const Point& v1, const VectorValue<Complex>& v2, const TensorValue<Complex> tensor) {

      Complex result = 0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          result += tensor(i, j) * v1(i) * v2(j);
        }
      }

      return result;
    }

    VectorValue<Complex> getVectorValue(Point p) {
      return VectorValue<Complex>(p);
    }
};

#endif
