# $Id$

topdir := ../..

include $(topdir)/Make.common


.PHONY: run

DIM := 1


run: diode.msh
	@echo "Running test02 ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad diode.in


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output tibercad


include $(topdir)/Make.rules
