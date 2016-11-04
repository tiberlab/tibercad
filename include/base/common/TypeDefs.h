// $Id$

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include "libmesh/id_types.h"
#include <limits.h>


//! A typedef to be used for numerical identifiers
/*!
 * We use the same type as the libMesh subdomain id type
 */
typedef libMesh::subdomain_id_type ID;

#ifndef INVALID_ID
# define INVALID_ID std::numeric_limits<ID>::max()
#endif



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
