# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

include $(topdir)/Make.common

target := modules

modules: libuptight module.conf
	@$(BUILDMODULE) module.conf

.PHONY: libuptight
libuptight:
	-@$(MAKE) --no-print-directory -C libuptight $(MAKECMDGOALS)
	@mkdir -p $(topdir)/$(ARCH)/modules/empirical_tb
	-@cp libuptight/libuptight.so $(topdir)/$(ARCH)/modules/empirical_tb



include $(topdir)/Make.rules
