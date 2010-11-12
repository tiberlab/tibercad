# $Id$

topdir := ../../../../..

MODULE_LDFLAGS = $(libdir)/modules/$(MODULE_NAME)/driftdiffusion$(shared_libext)

include $(topdir)/Make.module

$(foreach lib, $(srcfiles), $(eval $(call MODULE, $(lib))))


include $(topdir)/Make.rules
