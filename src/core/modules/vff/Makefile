# $Id: Makefile 2399 2011-02-23 15:43:57Z maufder $

topdir := ../../../..

MODULE_NAME := vff
export MODULE_NAME

include $(topdir)/Make.module


$(eval $(call MODULE, $(srcfiles)))


#target += models  #Models compilation temporary deactivated

.PHONY: models
models:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

clean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

distclean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)
	

include $(topdir)/Make.rules
