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

    enum Type {ElectricWall,  MagneticWall, AirWall};

  public:
    Type type;
    std::complex<double> alpha;

    bool isDirichle() const {
      return type == ElectricWall;
    }

    bool isMixed() const {
      return type == AirWall;
    }

  protected:
    MaxwellBoundaryProperties(const ModelOptions& options) : BoundaryProperties(options) {
      const std::string sType = options.get_option("type", "ElectricWall");
      if (sType == "ElectricWall") {
        type = ElectricWall;
      } else if (sType == "MagneticWall") {
        type = MagneticWall;
      } else {
        type = AirWall;
        alpha = std::complex<double>(options.get_option("em_alpha_r", 1.0), options.get_option("em_alpha_c", 1.0));
      }
    }
};

#endif
