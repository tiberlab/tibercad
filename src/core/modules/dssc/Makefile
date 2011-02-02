# $Id$

topdir := ../../../..

MODULE_NAME := dssc

include $(topdir)/Make.module


$(eval $(call MODULE, DSSC.C DSSCContact.C DSSCModel.C))
$(eval $(call MODULE, DSSCContact.C ))
$(eval $(call MODULE, DSSCModel.C))
$(eval $(call MODULE, DSSCGeneration.C))

MODULE_NAME := dssc_generation
$(eval $(call MODULE, DSSCGeneration.C))


include $(topdir)/Make.rules
