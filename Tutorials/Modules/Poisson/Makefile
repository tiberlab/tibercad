# $Id$

topdir := ../../..

MODULE_NAME := mypoisson
export MODULE_NAME

include $(topdir)/Make.module


srcfiles := MyPoisson.C MyPoissonModel.C PoissonBoundaryModel.C
models := DefaultModel.C PoissonDirichlet.C


$(eval $(call MODULE, $(srcfiles)))

$(foreach lib, $(models), $(eval $(call MODULE, $(lib))))


include $(topdir)/Make.rules
