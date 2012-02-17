#ifndef FUNCTION_INFO_H_
#define FUNCTION_INFO_H_

#include "ItemId.h"

class FunctionInfo {
  public:
    const static unsigned int MAX_ORDER = 10;

    enum ItemType {
      VERTEX, EDGE,
      SURFACE, VOLUME,
      UNDEFINED
    };

    ItemType itemType;
    unsigned int order;
    unsigned int itemId;

    ItemId globalItemId;
    bool isInterior;

    FunctionInfo() : itemType(UNDEFINED), order(0), itemId(0), isInterior(false) {
    }

    FunctionInfo(ItemType item, int o, int id) : itemType(item), order(o), itemId(id), isInterior(false) {
    }

    bool operator<(const FunctionInfo& foo) const {
      if (itemType < foo.itemType) {
        return true;
      } else if (itemType == foo.itemType) {
        if (order < foo.order) {
          return true;
        } else if (order == foo.order) {
          return globalItemId < foo.globalItemId;
        }
      }

      return false;
    }

    void setGlobalId(const ItemId& newGlobalId) {
      globalItemId = newGlobalId;
    }

    void setIsInteriour(bool flag) {
      isInterior = flag;
    }
};

#endif
