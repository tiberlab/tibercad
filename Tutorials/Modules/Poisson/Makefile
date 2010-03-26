# $Id$

topdir := ../../..

# subdirectories
MODULE_LDFLAGS := -L$(topdir)/lib -ltibercad

include $(topdir)/Make.module

MODULENAME=mypoisson

mainsrc := MyPoisson.C
srcfiles := MyPoissonModel.C PoissonBoundaryModel.C
models := DefaultModel.C
bdmodels := PoissonDirichlet.C


$(eval $(call MODULE, $(mainsrc), $(srcfiles)))
$(foreach lib, $(models), $(eval $(call SUBMODEL, $(lib))))
$(foreach lib, $(bdmodels), $(eval $(call SUBMODEL, $(lib))))


include $(topdir)/Make.rules
