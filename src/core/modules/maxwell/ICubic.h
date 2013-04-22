#ifndef ICUBIC_H_
#define ICUBIC_H_

#include "elem.h"
#include "point.h"
#include "MaxwellEquationsCommon.h"
#include "CubicEigenSystem.h"
#include "Constants.h"

// TODO
// How to name this class?

class ICubic {
  protected:
    Complex W1;
    double scaling;

  public:
    virtual bool addCData(CubicEigenSystem& system) {
    }

    virtual void init(CubicEigenSystem& system) {
      scaling = system.getEdgeDofMap(false)->getGeometryEx()->getScaling().get_length_scaling() *
          system.simulationInterface->get_environment().get_device().get_mesh_units() / Constants::c;
      Complex W1_scaled = W1 * scaling;

      CubicEigenSystem& cubicSystem = dynamic_cast<CubicEigenSystem&>(system);
      if (std::abs(W1_scaled - cubicSystem.getLambda1()) > 0.01 * std::abs(cubicSystem.getLambda1()) && std::abs(cubicSystem.getLambda1()) > 1e-10) {
        throw ModelErrorException("ICubic: different values of lambda1 are not supported!");
      }
      cubicSystem.setLambda1(W1_scaled);
    }
};

#endif
