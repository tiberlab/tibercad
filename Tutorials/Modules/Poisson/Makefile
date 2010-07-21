# $Id$

topdir := ../../..

MODULE_NAME := mypoisson
export MODULE_NAME

include $(topdir)/Make.module

$(eval $(call MODULE, $(srcfiles)))


target += models


.PHONY: models
models:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

clean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

distclean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)
	

include $(topdir)/Make.rules
