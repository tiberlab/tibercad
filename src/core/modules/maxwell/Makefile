# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

MODULE_NAME := maxwell

include $(topdir)/Make.module

$(eval $(call MODULE, MaxwellEquations.C CubicEigenSystem.C EdgeDofMap.C EigenSystem.C ScalarFEBase2D.C ScalarFEBase3D.C VectorFEBase2D.C VectorFEBase3D.C))

#MODULE_LIB =  $(module_dir)/$(MODULE_NAME)/maxwell$(shared_libext)

MODULE_NAME := maxwell_boundary

$(eval $(call MODULE, MaxwellBoundaryEquations.C, $(MODULE_LIB)))
#$(eval $(call MODULE, MaxwellBoundaryEquations.C))

include $(topdir)/Make.rules


