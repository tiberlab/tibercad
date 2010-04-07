# $Id$

topdir := ../../../..

MODULENAME := driftdiffusion
export MODULENAME

MODULE_LDFLAGS := -L$(topdir)/lib -ltibercad

include $(topdir)/Make.module


mainsrc := DriftDiffusion.C
srcfiles := $(filter-out $(mainsrc), $(srcfiles))


$(eval $(call MODULE, $(mainsrc), $(srcfiles)))


target += models

.PHONY: models
models:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

clean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)

distclean:
	@$(MAKE) --no-print-directory -C models $(MAKECMDGOALS)
	

include $(topdir)/Make.rules
