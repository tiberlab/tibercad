#ifndef EDGEDOFMAP_H_
#define EDGEDOFMAP_H_

#include <dof_map.h>
#include <mesh_base.h>
#include <map>
#include <set>
#include "IGeometryEx.h"
#include "FunctionInfo.h"
#include "VariableType.h"
#include "ItemId.h"

class EdgeDofMap {//TODO rename
  public:
    EdgeDofMap(IGeometryEx* geometryEx, std::vector<VariableType>& vars);//TODO &&

    void dof_indices(const Elem* const elem, std::vector<unsigned int>& di);

    void dof_indices(const Elem* const elem, std::vector<unsigned int>& di, unsigned int vn, bool resizeToZero = true);

    unsigned int getFreedomDegree();

    // Function info -> function index
    // For all variables
    // Interior functions are not stored
    std::vector<std::map<FunctionInfo, unsigned int>> reverseFunctionsIndeces;

    // return plevel for specified element and variable
    unsigned int getPOrder(const Elem* elem, unsigned int vn);

    IGeometryEx* getGeometryEx() {
      return geometryEx;
    }

    void check_directions(unsigned int varNum);

    // Element -> order
    // Stored for all elements
    std::vector<std::map<ItemId, int>> POrders;

  private:
    static const int CHECK_POINTS = 10;

    // Indeces for edge functions. function index -> Function info
    // Stored for all variables
    std::vector<std::map<unsigned int, FunctionInfo>> functionsIndeces;

    // For each variable store:
    //     For each element store all indeces
    // We use centroid of element as id
    std::vector<std::map<ItemId, std::vector<unsigned int>>> elemIndeces;

    unsigned int freedomDegree;

    std::vector<VariableType>& variables;

    IGeometryEx* geometryEx;

    std::vector<Point> constructPoints(Point p1, Point p2, int direction);
};

#endif /* EDGEDOFMAP_H_ */

