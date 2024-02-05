# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

include $(topdir)/Make.common

target := modules
#target := 

modules: libnegf module.conf
	@$(BUILDMODULE) module.conf >> test/double_barrier/Make_log_3.out 2>&1

.PHONY: libnegf clean
libnegf:
	-@$(MAKE) --no-print-directory -C libnegf $(MAKECMDGOALS) >> test/double_barrier/Make_log_2.out 2>&1
	@mkdir -p $(topdir)/$(ARCH)/modules/negf
	-@cp libnegf/libnegf.so $(topdir)/$(ARCH)/modules/negf

clean:
	cd libnegf; $(MAKE) clean


include $(topdir)/Make.rules
