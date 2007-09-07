# $Id$

topdir := ../..

include $(topdir)/Make.common


.PHONY: run

DIM := 1


run: bulk.msh
	@echo "Running test01 ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad bulk.in


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output tibercad


include $(topdir)/Make.rules
