/*
 * DofEdgeMap.cpp
 *
 *  Created on: Feb 10, 2011
 *      Author: paveryan
 */

#include "EdgeDofMap.h"
#include "VectorFEBase2D.h"
#include "VectorFunction.h"
#include "ElementUtils.h"
#include "FunctionInfo.h"
#include "IFEBase.h"
#include "ScalarFEBase2D.h"

using namespace std;

EdgeDofMap::EdgeDofMap(IGeometryEx* geometry, std::vector<VariableType>& vars) : variables (vars), geometryEx(geometry) {


  //End init variables.
  MeshBase& mesh = geometryEx->getMesh();
  mesh.find_neighbors();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  freedomDegree = 0;
  reverseFunctionsIndeces.resize(variables.size());
  functionsIndeces.resize(variables.size());

  for (unsigned int varNum = 0; varNum < variables.size(); varNum++) {
    std::map<ItemId, std::vector<unsigned int>> oneElemIndeces;

    IIFEBase* fe = variables[varNum].getFEbase();

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    for (; el != end_el; ++el) {
      const Elem* elem = *el;

      vector<unsigned int> elemIndecesVector;

      unsigned int order = getPOrder(elem, varNum);

      fe->reinit(elem, order);

      vector<FunctionInfo>& functionsInfo = fe->getFunctionsInfo();

      for (unsigned int j = 0; j < functionsInfo.size(); j++) {
        FunctionInfo info = functionsInfo[j];

        //std::cout << "Func \n"; flush(std::cout);
        //std::cout << "1\n"; flush(std::cout);
        geometry->setGlobalId(info, elem);
        bool excludeFunction = geometry->excludeFunction(elem, info);
        //std::cout << "excludeFunction: " << excludeFunction << "\n"; flush(std::cout);
        //std::cout << "info.isInterior: " << info.isInterior << "\n\n"; flush(std::cout);
        //std::cout << "3\n"; flush(std::cout);

        unsigned int functionIndex;
        if (info.isInterior) {
          // interior function
          functionIndex = freedomDegree;
          freedomDegree++;
        } else if (excludeFunction) {
          functionIndex = ElementUtils::INVALID_FUNCTION_ID;
        } else if (reverseFunctionsIndeces[varNum].find(info) == reverseFunctionsIndeces[varNum].end()) {
          //std::cout << "New GlobalId: " << info.globalItemId.id << "\n"; flush(std::cout);
          functionIndex = freedomDegree;
          freedomDegree++;

          reverseFunctionsIndeces[varNum].insert(make_pair(info, functionIndex));
          functionsIndeces[varNum].insert(make_pair(functionIndex, info));
        } else {
          functionIndex = reverseFunctionsIndeces[varNum][info];
        }

/*
        bool sourceFunction = geometry->sourceFunction(elem, info);
        if (sourceFunction && info.order == 0) {
          sourceIds.insert(functionIndex);
        }
*/
        elemIndecesVector.push_back(functionIndex);
      }

      oneElemIndeces.insert(make_pair(ItemId::get(elem), elemIndecesVector));
    }
    elemIndeces.push_back(oneElemIndeces);
    delete(fe);
  }
}

void EdgeDofMap::dof_indices(const Elem* const elem, vector<unsigned int>& di) {
  di.resize(0);
  for (unsigned int varNum = 0; varNum < variables.size(); varNum++) {
    dof_indices(elem, di, varNum, false);
  }
}

void EdgeDofMap::dof_indices(const Elem* const elem, vector<unsigned int>& di, unsigned int varNum, bool resizeToZero) {
  if (resizeToZero) {
    di.resize(0);
  }

  vector<unsigned int> indices = elemIndeces[varNum][ItemId::get(elem)];
  for (int i = 0; i < indices.size(); i++) {
    di.push_back(indices[i]);
  }
}

unsigned int EdgeDofMap::getFreedomDegree() {
  return freedomDegree;
}

unsigned int EdgeDofMap::getPOrder(const Elem* elem, unsigned int varNum) {
  while (POrders.size() <= varNum) {
    std::map<ItemId, int> tmp;
    POrders.push_back(tmp);
  }

  const Elem* temp = elem;

  while (temp != NULL && POrders[varNum].find(ItemId::get(temp)) == POrders[varNum].end()) {
    temp = temp->parent();
  }

  int order = (temp == NULL) ? variables.at(varNum).order : POrders[varNum][ItemId::get(temp)];

  POrders[varNum].insert(std::make_pair(ItemId::get(elem), order));

  return order;
}

void EdgeDofMap::check_directions(unsigned int varNum) {
/*  MeshBase::const_element_iterator el = geometryEx->getMesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = geometryEx->getMesh().active_local_elements_end();

  IVectorFEBase* fe = NULL;
  IScalarFEBase* fe_sc = NULL;
  if (geometryEx->getMesh().mesh_dimension() == 2) {
    if (variables[varNum].isVector) {
      fe = dynamic_cast<VectorFEBase2D*>(variables[varNum].getFEbase());
    } else {
      fe_sc = dynamic_cast<ScalarFEBase2D*>(variables[varNum].getFEbase());
    }
  } else if (geometryEx->getMesh().mesh_dimension() == 3) {
    if (variables[varNum].isVector) {
      fe = dynamic_cast<VectorFEBase3D*>(variables[varNum].getFEbase());
    } else {
      fe_sc = dynamic_cast<ScalarFEBase3D*>(variables[varNum].getFEbase());
    }
  } else {
    return;//do not need to check directions in 1D
  }

  std::map<std::pair<int, int>, std::vector<Point> > functionValues;
  std::map<std::pair<int, int>, std::vector<double> > functionValues_sc;

  for (; el != end_el; ++el)
  {
    const Elem* elem = *el;
    std::vector<unsigned int> dof_indices;
    EdgeDofMap::dof_indices (elem, dof_indices, varNum);

    std::vector<int> verticesIds;
    ElementUtils::getVertecesIds(elem, verticesIds);

    std::vector<std::vector<Point>> midPoints;
    if (verticesIds.size() == 4) {
      midPoints.push_back(constructPoints(Point(-1, -1), Point(1, -1), ElementUtils::getDirection(elem, 0)));
      midPoints.push_back(constructPoints(Point(1, -1), Point(1, 1), ElementUtils::getDirection(elem, 1)));
      midPoints.push_back(constructPoints(Point(1, 1), Point(-1, 1), ElementUtils::getDirection(elem, 2)));
      midPoints.push_back(constructPoints(Point(-1, 1), Point(-1, -1), ElementUtils::getDirection(elem, 3)));
    } else if (verticesIds.size() == 3) {
      midPoints.push_back(constructPoints(Point(0, 0), Point(1, 0), ElementUtils::getDirection(elem, 0)));
      midPoints.push_back(constructPoints(Point(1, 0), Point(0, 1), ElementUtils::getDirection(elem, 1)));
      midPoints.push_back(constructPoints(Point(0, 1), Point(0, 0), ElementUtils::getDirection(elem, 2)));
    }

    for (int i = 0; i < midPoints.size(); i++) {

      if (variables[varNum].isVector) {
        fe->reinit (elem, getPOrder(elem, varNum), &midPoints[i]);

        const std::vector<VectorFunction >& edge_phi = fe->getFunctions();

        int i1 = ElementUtils::getNext(i, verticesIds.size());
        ItemId edgeId = ItemId::get(elem->build_edge(verticesIds[i]));

        Point tangent = (*elem).point(i1) - (*elem).point(i);


        for (int j = 0; j < dof_indices.size(); j++) {// Iterate over functions
          std::vector<Point> newValues;
          for (int uu = 0; uu < 10; uu++) {
            newValues.push_back(edge_phi[j].phi[uu]);
          }

          if (functionValues.find(std::make_pair(edgeIndex, dof_indices[j])) == functionValues.end()) {
            functionValues[std::make_pair(edgeIndex, dof_indices[j])] = newValues;
          } else {
            std::vector<Point> oldValues = functionValues[std::make_pair(edgeIndex, dof_indices[j])];

            for (int uu = 0; uu < 10; uu++) {
              if (abs((oldValues[uu] - newValues[uu]) * tangent) > 1e-3 && dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
                std::cout << "Vfunction value is wrong: function " << dof_indices[j] << " on edge " << edgeIndex << " tangent: " << tangent(0) << " " << tangent(1) << "\n";
                std::cout << "                                  value1: " << oldValues[uu](0) << " " << oldValues[uu](1) <<
                                                          " value2: " << newValues[uu](0) << " " << newValues[uu](1) << " on point: " << uu << "\n";
                flush(std::cout);
                break;
              }
            }
          }
        }
      } else {
        fe_sc->reinit (elem, getPOrder(elem, varNum), &midPoints[i]);

        const std::vector<ScalarFunction >& phi = fe_sc->getFunctions();

        int i1 = ElementUtils::getNext(i, verticesIds.size());
        int edgeIndex = geometryEx.getEdgeIndex(verticesIds[i], verticesIds[i1]);

        for (int j = 0; j < dof_indices.size(); j++) {// Iterate over functions
          std::vector<double> newValues;
          for (int uu = 0; uu < 10; uu++) {
            newValues.push_back(phi[j].phi[uu]);
          }

          if (functionValues_sc.find(std::make_pair(edgeIndex, dof_indices[j])) == functionValues_sc.end()) {
            functionValues_sc[std::make_pair(edgeIndex, dof_indices[j])] = newValues;
          } else {
            std::vector<double> oldValues = functionValues_sc[std::make_pair(edgeIndex, dof_indices[j])];

            for (int uu = 0; uu < 10; uu++) {
              if (abs((oldValues[uu] - newValues[uu])) > 1e-3 && dof_indices[j] != ElementUtils::INVALID_FUNCTION_ID) {
                std::cout << " SFunction value is wrong: function " << dof_indices[j] << " on edge " << edgeIndex << "\n";
                std::cout << "                                  value1: " << oldValues[uu] << " value2: " << newValues[uu] << " on point: " << uu << "\n";
                flush(std::cout);
                break;
              }
            }
          }
        }
      }
    }
  }

  delete fe;
  delete fe_sc;*/
}

std::vector<Point> EdgeDofMap::constructPoints(Point p1, Point p2, int direction) {
  std::vector<Point> result;
  for (int i = 0; i <= EdgeDofMap::CHECK_POINTS; i++) {
    if (direction == 1) {
      result.push_back(p1 + (i+1.0)/(EdgeDofMap::CHECK_POINTS+2.0)*(p2 - p1)); // Some points between p1 and p2.
    } else {
      result.push_back(p2 + (i+1.0)/(EdgeDofMap::CHECK_POINTS+2.0)*(p1 - p2)); // Some points between p1 and p2.
    }
  }
  return result;
}
