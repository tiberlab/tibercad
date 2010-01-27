// $Id$

#ifndef _IDSET_H_
#define _IDSET_H_

#include "TypeDefs.h"
#include "HashSet.h"
#include <set>

//! A convenient typedef for a set of IDs
typedef std::set<ID> IDSet;

//! A convenient typedef for a hash set of IDs
typedef TiberCad::HashSet<ID>::Type IDHashSet;


#endif // _IDSET_H_
