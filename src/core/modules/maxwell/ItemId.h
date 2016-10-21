#ifndef ITEMID_H_
#define ITEMID_H_

#include "elem.h"
#include "point.h"

class ItemId {
  public:
    libMesh::Point id;

    ItemId() {
      id = 0;
    }

    ItemId(const libMesh::Point& iden) {
      id = iden;
    }

  public:
    bool operator<(const ItemId& foo) const {
      return id < foo.id;
    }

    bool operator==(const ItemId& foo) const {
      return id == foo.id;
    }

    static ItemId get(const libMesh::Elem* elem) {
      return ItemId(elem->centroid());
    }

    static ItemId get(const libMesh::Elem* elem, unsigned int vertex) {
      return ItemId(elem->point(vertex));
    }

    static ItemId get(const libMesh::Node* node) {
      return ItemId(*node);
    }
};

#endif
