/*
 * IGeometryEx.h
 *
 *  Created on: May 25, 2011
 *      Author: paveryan
 */

#ifndef IGEOMETRYEX_H_
#define IGEOMETRYEX_H_

#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Boundary.h"
#include "ElementSide.h"
#include "MaxwellBoundaryProperties.h"
#include "ElementUtils.h"
#include <mesh_base.h>
#include "Scaling.h"
#include "limits.h"
#include "FunctionInfo.h"
#include "Database.h"
#include "Material.h"
#include "PML.h"

class IGeometryEx {
  public:
    PML pml;
  protected:
    Scaling& scaling;
    static const double SCALE_LENGTH = 10;
    SimulationInterface* simulationInterface;
    std::set<ItemId> boundaryIds;

    //'Max' and 'Min' points of our device. (PML layers are included)
    Point minPoint;
    Point maxPoint;

  public:
    IGeometryEx(SimulationInterface* interface, Scaling& sc) : simulationInterface(interface), scaling(sc) {
      MeshBase& mesh = simulationInterface->get_mesh();

      {
        minPoint = Point(INT_MAX, INT_MAX);
        maxPoint = Point(INT_MIN, INT_MIN);

        const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        for (; el != end_el; ++el) {
          const Elem* elem = *el;

          for (unsigned int i = 0; i < elem->n_nodes(); i++) {
            const Node node = *(elem->get_node(i));

            minPoint(0) = std::min(minPoint(0), node(0));
            minPoint(1) = std::min(minPoint(1), node(1));
            maxPoint(0) = std::max(maxPoint(0), node(0));
            maxPoint(1) = std::max(maxPoint(1), node(1));
          }
        }

        //To avoid looking for this stupid error...
        if (mesh.mesh_dimension() == 1 && (maxPoint(1) - minPoint(1)) > 1e-5) {
          std::cout << "********************************************************************** \n";
          std::cout << "WARNING! Y_COORDINATE DIFFERS! PROBABLE YOU ARE TRYING CALULATE 1D WITH Y-MESH!\n";
        }

        double structureDiameter = (maxPoint - minPoint).size();
        std::cout << "Structure diameter: " << structureDiameter << "\n";
        scaling.set_length_scaling(structureDiameter / SCALE_LENGTH);
        std::cout << "Using scale: " << scaling.get_length_scaling() << "\n";
        flush(std::cout);
      }


      mesh.find_neighbors();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
      MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
      for (; el != end_el; ++el) {
        const Elem* elem = *el;

        for (unsigned int i = 0; i < elem->n_sides(); i++) {
          if (excludeSideFunction(elem, i)) {
            const AutoPtr<Elem> sideElem = elem->build_side(i);

            //std::cout << "Excluded " << i << " " << sideElem->centroid() << "\n";

            boundaryIds.insert(ItemId::get(sideElem.get()));
            for (int j = 0; j < sideElem->n_nodes(); j++) {
              boundaryIds.insert(ItemId::get(sideElem->get_node(j)));
            }
            if (mesh.mesh_dimension() == 3) {
              for (int j = 0; j < sideElem->n_edges(); j++) {
                boundaryIds.insert(ItemId::get(sideElem->build_edge(j).get()));
              }
            }
          }
        }
      }

      pml.init(interface);
      pml.allMaxPoint = maxPoint;
      pml.allMinPoint = minPoint;
      //std::cout << "BE " << boundaryIdsEdges.size() << "\n"; flush(std::cout);
    }

    virtual ~IGeometryEx() {}

    virtual Scaling& getScaling() {
      return scaling;
    }

    virtual MeshBase& getMesh() {
      return simulationInterface->get_mesh();
    }

    virtual MaxwellBoundaryProperties* getBoundaryProperties(const Elem* elem, unsigned int side) {
      if (!isSideBoundary(elem, side)) {
        return NULL;
      }

      Boundary* bd = simulationInterface->get_environment().get_boundary(ElementSide(elem,side));

      if (bd == NULL || (bd->get_boundary_properties( simulationInterface->get_id() ) == NULL )) {
        return NULL;
      }

      return dynamic_cast<MaxwellBoundaryProperties*>(bd->get_boundary_properties(simulationInterface->get_id()));
    }

    virtual bool excludeSideFunction(const Elem* elem, unsigned int side) {
      if (!isSideBoundary(elem, side)) {
        return false;
      }

      MaxwellBoundaryProperties* properties = getBoundaryProperties(elem, side);

/*
      std::cout << "Got props " <<  properties << "\n"; flush(std::cout);
      if (properties != NULL) {
        std::cout << "Got props " <<  properties->isDirichle() << "\n"; flush(std::cout);
      }
*/

      return properties == NULL || properties->isDirichle();
    }

    virtual bool excludeFunction(const Elem* elem, FunctionInfo& info) {
      return boundaryIds.count(info.globalItemId) != 0;
    }

    virtual bool isSideBoundary(const Elem* elem, unsigned int side) {
      return elem->neighbor(side) == NULL;
    }

    //

    virtual void setGlobalId(FunctionInfo& info, const Elem* elem) {
      if (getMesh().mesh_dimension() == 3) {
        info.setIsInteriour(info.itemType == FunctionInfo::VOLUME);
      } else if (getMesh().mesh_dimension() == 2) {
        info.setIsInteriour(info.itemType == FunctionInfo::SURFACE);
      } else {
        info.setIsInteriour(info.itemType == FunctionInfo::EDGE);
      }

      if (!info.isInterior) {
        ItemId globalItemId;
        if (info.itemType == FunctionInfo::VERTEX) {
          globalItemId = ItemId::get(elem, info.itemId);
        } else if (info.itemType == FunctionInfo::EDGE) {
          globalItemId = ItemId::get(elem->build_edge(info.itemId).get());
        } else if (info.itemType == FunctionInfo::SURFACE) {
          globalItemId = ItemId::get(elem->build_side(info.itemId).get());
        }
        info.setGlobalId(globalItemId);
      }
    }
};

#endif /* IGEOMETRYEX_H_ */
