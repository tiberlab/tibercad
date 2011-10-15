# $Id$

topdir := ../../../..


all: module.conf
	@export TIBERCADROOT=$(topdir) ;  $(topdir)/i686-linux/bin/build_module module.conf


.PHONY: models
models:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

clean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

distclean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)
	

