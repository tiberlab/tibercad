#ifndef _ELEMENT_UTILS_
#define _ELEMENT_UTILS_

#include <vector>
#include "point.h"
#include <TypeDefs.h>
#include "elem.h"

class ElementUtils
{
  public:
    static const unsigned int INVALID_FUNCTION_ID = INVALID_ID;

    enum ElementType {
      TRI, QUAD
    };

    static ElementType getType(const Elem* elem) {
      return TRI;
    }

    static void getVertecesIds(const Elem* elem, std::vector<int>& result) {
      int nodesCount = elem->n_nodes();

      for (int i = 0; i < nodesCount; i++) {
        if (elem->is_vertex(i)) {
          result.push_back(elem->get_node(i)->id());
        }
      }
    }

    static unsigned int getNext(const unsigned int i, const unsigned int size) {
      return (i  == size - 1) ? 0 : (i + 1);
    }

    static unsigned int getPrevious(const unsigned int i, const unsigned int size) {
      return i == 0 ? (size - 1) : (i - 1);
    }

    /* Rewrite the following methods if you want support i.e. TRI6 or QUAD9 elements.
     * Thow it will work slower it will work...*/

    static int getDirection(const Elem* elem, const unsigned int index) {
      int i1 = (index == elem->n_nodes() - 1) ? 0 : (index + 1);
      return ( elem->get_node(index)->id() < elem->get_node(i1)->id() ) ? 1 : -1;
    }

    static int getDirection(const Elem* elem, const unsigned int index1, const unsigned index2) {
      return ( elem->get_node(index1)->id() < elem->get_node(index2)->id() ) ? 1 : -1;
    }

    static unsigned int getVertecesCount(const Elem* elem) {
      return elem->n_nodes();
    }

    static unsigned int getVertexId(const Elem* elem, const unsigned int i) {
      return elem->get_node(i)->id();
    }

    static std::vector<std::pair<unsigned int, unsigned int>> getSortedEdge(const Elem* elem, unsigned int edgeNum) {
      std::vector<std::pair<unsigned int, unsigned int>> result;

      for (unsigned int i = 0; i < elem->n_nodes(); i++) {
        if (elem->is_node_on_edge(i, edgeNum)) {
          result.push_back(std::make_pair(i, elem->get_node(i)->id()));
        }
      }

      std::sort(result.begin(), result.end(), comparePair);
      return result;
    }

    /**
     * Result is map<vertex element index, global vertex index>, where vertices are sort according to global indices.
     */
    static std::vector<std::pair<unsigned int, unsigned int>> getSortedSide(const Elem* elem, unsigned int sideNum) {
      std::vector<std::pair<unsigned int, unsigned int>> result;

      for (unsigned int i = 0; i < elem->n_nodes(); i++) {
        if (elem->is_node_on_side(i, sideNum)) {
          result.push_back(std::make_pair(i, elem->get_node(i)->id()));
        }
      }

      std::sort(result.begin(), result.end(), comparePair);
      return result;
    }

    static bool comparePair(std::pair<unsigned int, unsigned int> pair1, std::pair<unsigned int, unsigned int> pair2) {
      return pair1.second < pair2.second;
    }
};

#endif
