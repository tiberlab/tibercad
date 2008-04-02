# $Id: Makefile 473 2007-09-07 15:06:21Z maufder $

topdir := ../..

inputfile = test_03.tib


.PHONY: run clean distclean


run: 
	@echo "Running test in $(notdir $(PWD)) ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad $(inputfile)


clean:
	@rm -rf *~ *.xda *.xta *.ucd tibercad output


