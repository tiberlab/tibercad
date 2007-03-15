#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <utility>
#include <map>
#include <vector>

class Elem;


//! An element side
typedef std::pair<const Elem*, unsigned int> ElementSide;


//! A typedef to be used for numerical identifiers
typedef unsigned int ID;


//! A map which contains all nodes belonging to boundary conditions
typedef std::map<ID, std::vector<ID> > BoundaryNodeMap;


//! To ignore unused variables
/*!
 * If you have a variable \a which is unused you can call
 * \c ignore_unused_variable(a) so that the compiler gives no warning
 */
template <typename T>
inline
void ignore_unused_variable(T)
{
}

#endif // _TYPEDEFS_H_
