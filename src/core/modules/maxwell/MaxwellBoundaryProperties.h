// $Id: MaxwellEquations.h 1793 2010-02-10 08:51:09Z maufder $

#ifndef _MAXWELLBOUNDARYPROPERTIES_H_
#define _MAXWELLBOUNDARYPROPERTIES_H_

#include "BoundaryProperties.h"
#include <complex>

//! Class to solve Maxwell equations
class MaxwellBoundaryProperties : public BoundaryProperties
{
  public:
    static MaxwellBoundaryProperties* create(const ModelOptions& options) {
      return new MaxwellBoundaryProperties(options);
    }

    enum Type {ElectricWall,  MagneticWall, Source};

  public:
    Type type;
    int direction;
    double power;

    bool isDirichle() const {
      return type == ElectricWall;
    }

    bool isSource() const {
      return type == Source;
    }

  protected:
    MaxwellBoundaryProperties(const ModelOptions& options) : BoundaryProperties(options) {
      const std::string sType = options.get_option("type", "ElectricWall");
      if (sType == "ElectricWall") {
        type = ElectricWall;
      } else if (sType == "MagneticWall") {
        type = MagneticWall;
      } else {
        type = Source;

        direction = options.get_option("direction", -1);
        power = options.get_option("power", 1);
      }
    }
};

#endif
