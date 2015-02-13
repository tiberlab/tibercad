# $Id: Makefile 3777 2014-02-07 11:21:14Z maufder $

topdir := ../..

inputfile = quantum_dot_atom.tib


.PHONY: run clean distclean


run: 
	@echo "Running test in $(notdir $(PWD)) ..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@./tibercad $(inputfile)

#@if [ -e output/eigv.dat ] ; then \
#  echo "# n vb1 vb2 cb1 cb2"|cat - output/eigv.dat > output/eigv_res.dat; \
#fi





clean:
	@rm -rf *~ *.xda *.xta *.ucd tibercad output


