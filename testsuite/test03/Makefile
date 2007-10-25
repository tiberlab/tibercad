# $Id: Makefile 473 2007-09-07 15:06:21Z maufder $

topdir := ../..

include $(topdir)/Make.common


.PHONY: run

DIM := 1


run: diode.msh
	@echo "Running test03 ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad diode.in


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh
	@rm -rf output tibercad


include $(topdir)/Make.rules
