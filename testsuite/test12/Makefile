# $Id$

topdir := ../..

inputfile = bulk.tib


.PHONY: run clean distclean


run: 
	@echo "Running test in $(notdir $(PWD)) ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad $(inputfile)


clean:
	@rm -rf *~  tibercad output


