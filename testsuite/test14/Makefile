# $Id: Makefile 564 2007-11-23 09:02:50Z maufder $

topdir := ../..

inputfile = put_your_filename_here.tib


.PHONY: run clean distclean


run: 
	@echo "Running test in $(notdir $(PWD)) ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad $(inputfile)


clean:
	@rm -rf *~ *.xda *.xta *.ucd tibercad output


