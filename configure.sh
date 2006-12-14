#!/bin/sh

export CXX=g++-4.0.0

./configure --with-boost-prefix=/usr/pack/boost-1.33.1-ma \
  --with-libmesh-prefix=/usr/pack/libmesh-0.5.0_20051010cvs-ma
