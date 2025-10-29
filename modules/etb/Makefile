# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

include $(topdir)/Make.common

target := modules

modules: libuptight module.conf
	@$(BUILDMODULE) module.conf

.PHONY: libuptight distclean

libuptight:
	-@$(MAKE) --no-print-directory -C libuptight $(MAKECMDGOALS)
	@mkdir -p $(topdir)/$(ARCH)/modules/empirical_tb
	-@cp libuptight/libuptight.so $(topdir)/$(ARCH)/modules/empirical_tb

distclean:
	cd libuptight/src/api; $(MAKE) clean
	cd libuptight; $(MAKE) distclean
	$(MAKE) clean


include $(topdir)/Make.rules
