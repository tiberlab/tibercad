#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <utility>
#include <map>

class Elem;

typedef std::pair<const Elem*, unsigned int> ElemSide;

typedef std::map<const std::string, std::string> ModelOptions;

typedef unsigned int ID;


#endif // _TYPEDEFS_H_
