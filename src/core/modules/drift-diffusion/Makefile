# $Id$

topdir := ../../../..

# subdirectories
subdirs		:= models


include $(topdir)/Make.module


mainsrc := DriftDiffusion.C
srcfiles := $(filter-out $(mainsrc), $(srcfiles))


$(eval $(call MODULE, $(mainsrc), $(srcfiles)))


include $(topdir)/Make.rules
