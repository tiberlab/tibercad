#!/usr/bin/python

configure_options = [
  '--known-level1-dcache-size=32768',
  '--known-level1-dcache-linesize=32',
  '--known-level1-dcache-assoc=2',
  '--known-memcmp-ok=1',
 '--known-endian=little',
  '--known-sizeof-char=1',
  '--known-sizeof-void-p=4',
  '--known-sizeof-short=2',
  '--known-sizeof-int=4',
  '--known-sizeof-long=4',
  '--known-sizeof-long-long=8',
  '--known-sizeof-float=4',
  '--known-sizeof-double=8',
  '--known-sizeof-size_t=4',
  '--known-bits-per-byte=8',
  '--with-matlab=0',
  '--with-shared-libraries=0',
  '--with-dynamic-loading=0',
  '--with-mpi=0',
  '--with-cc=i686-w64-mingw32-gcc',
  '--with-cxx=i686-w64-mingw32-g++',
  '--with-fc=i686-w64-mingw32-gfortran',
  '--with-ranlib=i686-w64-mingw32-ranlib',
  '--with-cpp=i686-w64-mingw32-cpp',
  '--with-shared-ld=i686-w64-mingw32-ld',
  '--with-ar=i686-w64-mingw32-ar',
  '--with-clanguage=C++',
  '--download-f-blas-lapack=1',
  '--with-batch=1',
  '--with-x=0',
  '--with-host-vendor=w64',
  '--with-host-os=mingw32',
  '--with-debugging=0',
  '--with-scalar-type=complex',
  '--CFLAGS=-g -fexceptions',
  '--FFLAGS=-fexceptions -fexceptions',
  '--CXXFLAGS=-g',
  '--LDFLAGS=',
]
if __name__ == '__main__':
  import os
  import sys
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
