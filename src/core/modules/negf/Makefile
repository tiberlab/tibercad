# $Id: Makefile 2309 2011-02-07 21:39:43Z maufder $

topdir := ../../../..

MODULE_NAME := negf

include $(topdir)/Make.module

.PHONY: libnegf
target += libnegf

MODULE_LDFLAGS += -L$(libdir) -lnegf

$(eval $(call MODULE, Negf.C ))


libnegf:       
	$(MAKE) -C libnegf




include $(topdir)/Make.rules
