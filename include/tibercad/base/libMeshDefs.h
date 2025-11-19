// $Id$

#ifndef _LIBMESHDEFS_H_
#define _LIBMESHDEFS_H_


//! Definitions for easier use of libMesh stuff

#define LIBMESHCLASS(name) namespace libMesh { class name; } using libMesh::name
#define USELIBMESHTYPE(name) using libMesh::name

// these are needed almost everywhere
LIBMESHCLASS(Elem);
LIBMESHCLASS(Point);
LIBMESHCLASS(Node);
LIBMESHCLASS(MeshBase);

// they have to correspond to the libMesh typedefs!
namespace libMesh
{
typedef double Real;
typedef Real Number;
}

USELIBMESHTYPE(Real);
USELIBMESHTYPE(Number);

#endif // _LIBMESHDEFS_H_
