# $Id: Makefile 2212 2010-12-07 18:53:51Z gromano $

topdir := ../../../..

MODULE_NAME := molecular-dynamics
export MODULE_NAME


include $(topdir)/Make.module

$(eval $(call MODULE,$(srcfiles),$(wildcard *.a)))

#target += models


#.PHONY: models
models:
	@$(MAKE) --no-print-directory -C $(MAKECMDGOALS)

clean:
	@$(MAKE) --no-print-directory -C $(MAKECMDGOALS)

distclean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)
	

include $(topdir)/Make.rules
