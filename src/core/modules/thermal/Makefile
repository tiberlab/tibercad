# $Id: Makefile 3004 2011-10-17 22:03:35Z maufder $

topdir := ../../../..

include $(topdir)/Make.common

target := modules

modules: module.conf
	@$(topdir)/bin/build_module module.conf


clean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

distclean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)
	


include $(topdir)/Make.rules
