# $Id$

topdir := ../..

include $(topdir)/Make.common


.PHONY: run

DIM := 2

all: run


run:
	@echo "Running example in $(DIM)D..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad device.in


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output tibercad


include $(topdir)/Make.rules
