# $Id$

topdir := ../../../..

# subdirectories
#subdirs	:= models

MODULE_LDFLAGS := -L$(topdir)/lib -ltibercad

include $(topdir)/Make.module

MODULENAME=driftdiffusion

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
