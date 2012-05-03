# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

MODULE_NAME := maxwell

include $(topdir)/Make.module

$(eval $(call MODULE, MaxwellEquations.C CubicEigenSystem.C EdgeDofMap.C EigenSystem.C VectorSystem.C ScalarFEBase2D.C ScalarFEBase3D.C VectorFEBase2D.C VectorFEBase3D.C))
MODULE_NAME := maxwell_boundary

$(eval $(call MODULE, MaxwellBoundaryEquations.C EdgeDofMap.C VectorLinearSystem.C VectorSystem.C ScalarFEBase2D.C ScalarFEBase3D.C VectorFEBase2D.C VectorFEBase3D.C))
include $(topdir)/Make.rules


