# $Id$

topdir := ../../../../..

MODULE_LDFLAGS = $(topdir)/lib/modules/$(MODULE_NAME)/driftdiffusion.so

include $(topdir)/Make.module

$(foreach lib, $(srcfiles), $(eval $(call MODULE, $(lib))))


include $(topdir)/Make.rules
