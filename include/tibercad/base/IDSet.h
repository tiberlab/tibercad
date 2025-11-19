// $Id$

#ifndef _IDSET_H_
#define _IDSET_H_

#include "tibercad/base/TypeDefs.h"
#include "tibercad/base/HashSet.h"
#include <set>

//! A convenient typedef for a set of IDs
typedef std::set<ID> IDSet;

//! A convenient typedef for a hash set of IDs
typedef HashSet<ID>::Type IDHashSet;


#endif // _IDSET_H_
