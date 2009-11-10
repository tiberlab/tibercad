# $Id$

topdir := ../..

include $(topdir)/Make.common


.PHONY: run

DIM := 2

all: run


run: mosfet.msh
	@echo "Running example MOSFET..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad mosfet.tib


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output tibercad


include $(topdir)/Make.rules
