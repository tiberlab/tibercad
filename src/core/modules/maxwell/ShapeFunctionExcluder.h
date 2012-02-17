#ifndef _SHAPE_FUNCTION_EXCLUDER_
#define _SHAPE_FUNCTION_EXCLUDER_

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Boundary.h"
#include "ElementSide.h"
#include "MaxwellBoundaryProperties.h"

class ShapeFunctionExcluder {

  public:
    SimulationInterface* simulationInterface;

    bool excludeEdgeFunction(const Elem* elem, int side) {
      Boundary* bd = simulationInterface->get_environment().get_boundary(ElementSide(elem,side));

      if (bd == NULL || (bd->get_boundary_properties( simulationInterface->get_id() ) == NULL )) {
        return true; // By default if BM is not specifed >>> electric_wall
      }

      MaxwellBoundaryProperties* properties = dynamic_cast<MaxwellBoundaryProperties*>(bd->get_boundary_properties(simulationInterface->get_id()));

      return properties->isDirichle();
    }
};

#endif
