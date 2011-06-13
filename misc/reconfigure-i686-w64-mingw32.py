#!/usr/bin/python

configure_options = [
 '--with-memcmp-ok', '--with-endian=little',
  '--sizeof_char=1',
  '--sizeof_void_p=4',
  '--sizeof_short=2',
  '--sizeof_int=4',
  '--sizeof_long=4',
  '--sizeof_long_long=8',
  '--sizeof_float=4',
  '--sizeof_double=8',
  '--sizeof_size_t=4',
 '--bits_per_byte=8',
  '--configModules=PETSc.Configure', '--optionsModule=PETSc.compilerOptions', '--with-matlab=0', '--with-shared=0', '--with-dynamic=0', '--with-mpi=0', '--with-cc=i686-w64-mingw32-gcc', '--with-cxx=i686-w64-mingw32-g++', '--with-fc=i686-w64-mingw32-gfortran', '--with-ranlib=i686-w64-mingw32-ranlib', '--with-cpp=i686-w64-mingw32-cpp', '--with-clanguage=C++', '--with-batch=1', '--with-x=0', '--with-host-vendor=w64', '--with-host-os=mingw32', '--download-f-blas-lapack=1', '--with-debugging=0', '--CFLAGS= -fexceptions', '--FFLAGS= -fexceptions', '--CXXFLAGS=', '--LDFLAGS='
]
if __name__ == '__main__':
  import os
  import sys
  sys.path.insert(0, os.path.abspath(os.path.join('config')))
  import configure
  configure.petsc_configure(configure_options)
