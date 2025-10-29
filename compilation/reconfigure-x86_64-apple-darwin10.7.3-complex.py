#!/usr/bin/python

configure_options = [
 '--with-memcmp-ok',  '--sizeof_char=1',
  '--sizeof_void_p=8',
  '--sizeof_short=2',
  '--sizeof_int=4',
  '--sizeof_long=8',
  '--sizeof_long_long=8',
  '--sizeof_float=4',
  '--sizeof_double=8',
  '--sizeof_size_t=8',
 '--bits_per_byte=8',
  '--configModules=PETSc.Configure', '--optionsModule=PETSc.compilerOptions', '--with-matlab=0', '--with-shared=0', '--with-dynamic=0', '--with-mpi=0', '--with-cc=gcc', '--with-cxx=g++', '--with-fc=gfortran', '--with-ranlib=ranlib', '--with-cpp=cpp', '--with-shared-ld=ld', '--with-ar=ar', '--with-clanguage=C++', '--download-f-blas-lapack=1', '--with-batch=1', '--with-x=0', '--with-debugging=0', '--with-scalar-type=complex', '--CFLAGS=-fPIC  -fexceptions', '--FFLAGS=-fPIC  -fexceptions', '--CXXFLAGS=-fPIC ', '--LDFLAGS='
]
if __name__ == '__main__':
  import os
  import sys
  sys.path.insert(0, os.path.abspath(os.path.join('config')))
  import configure
  configure.petsc_configure(configure_options)
