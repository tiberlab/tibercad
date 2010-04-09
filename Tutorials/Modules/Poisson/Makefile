# $Id$

topdir := ../../..

MODULE_NAME := mypoisson
export MODULE_NAME

include $(topdir)/Make.module


mainsrc := MyPoisson.C
srcfiles := MyPoissonModel.C PoissonBoundaryModel.C
models := DefaultModel.C PoissonDirichlet.C


$(eval $(call MODULE, $(mainsrc), $(srcfiles)))

$(foreach lib, $(models), $(eval $(call MODULE, $(lib))))


include $(topdir)/Make.rules
