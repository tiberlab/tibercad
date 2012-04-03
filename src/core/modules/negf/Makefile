# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

include $(topdir)/Make.common

target := modules

modules: libnegf module.conf
	@$(BUILDMODULE) module.conf

.PHONY: libnegf
libnegf:
	-@$(MAKE) --no-print-directory -C libnegf $(MAKECMDGOALS)



include $(topdir)/Make.rules
