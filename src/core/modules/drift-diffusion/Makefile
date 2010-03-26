# $Id$

topdir := ../../../..

# subdirectories
subdirs	:= models
MODULE_LDFLAGS := -L$(topdir)/lib -ltibercad

include $(topdir)/Make.module

MODULENAME=driftdiffusion

mainsrc := DriftDiffusion.C
srcfiles := $(filter-out $(mainsrc), $(srcfiles))


$(eval $(call MODULE, $(mainsrc), $(srcfiles)))


include $(topdir)/Make.rules
