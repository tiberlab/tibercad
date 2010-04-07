# $Id$

topdir := ../../../../..

#MODULE_LDFLAGS = -L$(topdir)/lib/tibermodels -l:driftdiffusion.so -L$(topdir)/lib -ltibercad
MODULE_LDFLAGS = $(topdir)/lib/tibermodels/driftdiffusion.so -L$(topdir)/lib -ltibercad

include $(topdir)/Make.module

$(foreach lib, $(srcfiles), $(eval $(call SUBMODEL, $(lib))))


include $(topdir)/Make.rules
