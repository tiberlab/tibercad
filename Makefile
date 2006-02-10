# $Id$
#

topdir := .
dirs := src doc examples

include Make.common


.PHONY: all doc clean



all:
	@$(MAKE) -C src all
	

clean:
	@for i in $(dirs);	\
	do			\
	  $(MAKE) -C $$i clean;	\
	done


examples: all
	@$(MAKE) -C examples


doc:
	@$(MAKE) -C doc



