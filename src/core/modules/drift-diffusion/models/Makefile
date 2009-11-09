# $Id$

topdir := ../../../../..

MODULE_LDFLAGS = -L$(topdir)/lib/tibermodels -ldriftdiffusion -L$(topdir)/lib -ltibercad
include $(topdir)/Make.module


$(foreach lib, $(srcfiles), $(eval $(call SUBMODEL, $(lib))))


include $(topdir)/Make.rules
