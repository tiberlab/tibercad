# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

include $(topdir)/Make.module

MODULE_NAME := maxwell
$(eval $(call MODULE, MaxwellEquations.C CubicEigenSystem.C EdgeDofMap.C EigenSystem.C VectorSystem.C ScalarFEBase2D.C ScalarFEBase3D.C VectorFEBase2D.C VectorFEBase3D.C OpticPropsModel.C))

MODULE_NAME := maxwell_boundary
$(eval $(call MODULE, MaxwellBoundaryEquations.C EdgeDofMap.C VectorLinearSystem.C VectorSystem.C ScalarFEBase2D.C ScalarFEBase3D.C VectorFEBase2D.C VectorFEBase3D.C OpticPropsModel.C))


MODULE_LIB := $(module_dir)/$(MODULE_NAME)/maxwell_boundary$(shared_libext)

$(eval $(call MODULE, OpticPropsWDependant.C, $(MODULE_LIB)))

include $(topdir)/Make.rules
