# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../..

include $(topdir)/Make.common

target := modules
#target := 

modules: libnegf module.conf
	@$(BUILDMODULE) module.conf 

.PHONY: libnegf clean
libnegf:
	@$(MAKE) --no-print-directory -C libnegf
	@mkdir -p $(topdir)/$(ARCH)/modules/negf
	@cp -f libnegf/libnegf.so $(topdir)/$(ARCH)/modules/negf/

clean:
	@$(MAKE) --no-print-directory -C libnegf clean


include $(topdir)/Make.rules
