# build directory (relative)
BUILDDIR=`pwd`/../extern

# PETSc version
PETSCVERSION="2.3.3"
PATCHLEVEL="p13"

# libMesh version
LIBMESHVERSION=0.6.1

# whether to compile debug mode or not (y/n)
DEBUG=n

# architecture
ARCH=cygwin

OTHEROPTS="--with-boost-regex-lib=/usr/lib/libboost_regex-gcc-mt.a --with-boost-filesystem-lib=/usr/lib/libboost_filesystem-gcc-mt.a"
OTHEROPTS="${OTHEROPTS} --enable-license-check"
