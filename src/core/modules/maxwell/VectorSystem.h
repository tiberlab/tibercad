#ifndef __VECTOR_SYSTEM_H__
#define __VECTOR_SYSTEM_H__

#include "libmesh_config.h"
#include "equation_systems.h"
#include "system.h"

#include "EdgeDofMap.h"
#include "IGeometryEx.h"
#include "VariableType.h"

#include <vector>


using namespace libMeshEnums;


class VectorSystem {
  public:

    /**
     * Constructor.
     */
    VectorSystem ();

    /**
     * Destructor.
     */
    virtual ~VectorSystem ();

     virtual EdgeDofMap* getEdgeDofMap(const bool construct = false);


    unsigned int addVariable(const unsigned int order, const bool isVector, const unsigned int extraQOrder = 0, const bool inplane = false);

    const VariableType& getVariableType(const unsigned int i) const;

    unsigned int getVariablesCount() const;

    virtual void dof_indices(const Elem* const elem, std::vector<unsigned int>& di);

    virtual double getGeomScaling();


    static void countNonZeros(std::map<std::pair<int, int>, libMesh::Complex>& values, int* result, unsigned int size);

  protected:
    virtual unsigned int initSystemSize();

    /**
     * Initializes the member data fields associated with
     * the system, so that, e.g., \p assemble() may be used.
     */
    EdgeDofMap* edgeDofMap;

    IGeometryEx* geometryEx;

  public:
    //TODO Should not be here!
    SimulationInterface* simulationInterface;

    std::vector<VariableType> variables;
  
    IGeometryEx* getGeometryEx() {
      return geometryEx;
    }
};

#endif
