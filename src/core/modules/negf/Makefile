# $Id: Makefile 2309 2011-02-07 21:39:43Z maufder $

topdir := ../../../..

MODULE_NAME := negf
export MODULE_NAME

include $(topdir)/Make.module

$(eval $(call MODULE, $(srcfiles)))

target += libnegf
target += models

MODULE_LDFLAGS += -L$(libdir) -lnegf

.PHONY: models
models:
	$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)


.PHONY: libnegf
libnegf:
	$(MAKE) -C libnegf


include $(topdir)/Make.rules
