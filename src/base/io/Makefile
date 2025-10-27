# $Id$

topdir := ../../..

subdirs := Read_ISE_grid 

include $(topdir)/Make.common

csrcfiles := $(wildcard *.c)
cobjects := $(patsubst %.c, .$(ARCH)/%$(obj-suffix), $(csrcfiles))
target += $(cobjects)

ifeq ($(target_os),win32)
target += .$(ARCH)/dlfcn$(obj-suffix)

CPPFLAGS += -I$(topdir)/external/dlfcn-win32-r19

.$(ARCH)/dlfcn$(obj-suffix): $(topdir)/external/dlfcn-win32-r19/dlfcn.c
	@echo "Compiling win32 compatibility layer for <dlfcn.h> ..."
	@$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<
endif

include $(topdir)/Make.rules

#print a numbered list of variables in grace output
CPPFLAGS += -DNUM_VAR

