// Local includes
#include "VectorSystem.h"


// ------------------------------------------------------------

VectorSystem::VectorSystem () : edgeDofMap (NULL) {
}

VectorSystem::~VectorSystem () {
  if (edgeDofMap != NULL) {
    delete edgeDofMap;
    delete geometryEx;
  }
}

unsigned int VectorSystem::initSystemSize () {
  return getEdgeDofMap(true)->getFreedomDegree();
}

EdgeDofMap* VectorSystem::getEdgeDofMap(const bool construct) {
  if (construct || edgeDofMap == NULL) {
    if (edgeDofMap != NULL) {
      delete(edgeDofMap);
      delete geometryEx;
    }
    //TODO mb use init-reinit
    geometryEx = new IGeometryEx(simulationInterface, simulationInterface->get_scaling());

    //TODO bad:
    for (int i = 0; i < variables.size(); i++) {
      variables[i].scaling = simulationInterface->get_scaling().get_length_scaling();
      variables[i].dimension = simulationInterface->get_mesh().mesh_dimension();
    }

    edgeDofMap = new EdgeDofMap(geometryEx, variables);
  }
  return edgeDofMap;
}

unsigned int VectorSystem::addVariable(const unsigned int order, const bool isVector, const unsigned int extraQOrder, const bool inplane) {
  variables.push_back(VariableType(order, isVector, 0/* This will be modified later!*/, 1.0, extraQOrder, inplane));
  return variables.size() - 1;
}

const VariableType& VectorSystem::getVariableType(const unsigned int i) const {
  return variables[i];
}

unsigned int VectorSystem::getVariablesCount() const {
  return variables.size();
}

void VectorSystem::dof_indices(const Elem* const elem, std::vector<unsigned int>& di) {
  getEdgeDofMap(false)->dof_indices(elem, di);
}

void VectorSystem::countNonZeros(std::map<std::pair<int, int>, libMesh::Complex>& values, int* result, unsigned int size) {
  for (int i = 0; i < size; i++) {
    result[i] = 0;
  }

  std::map<std::pair<int, int>, libMesh::Complex>::iterator it = values.begin();
  for (; it != values.end(); it++) {
    std::pair<std::pair<int, int>, libMesh::Complex> pair = *it;
    result[pair.first.first]++;
  }
}

double VectorSystem::getGeomScaling() {
  return geometryEx->getScaling().get_length_scaling() * geometryEx->getScaling().get_length_scaling();
}
