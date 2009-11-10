# $Id$

topdir := ../../..

include $(topdir)/Make.common

targets := run1d run2d

infile := diode.in

.PHONY: run

DIM := 2

all: run

run1d: DIM := 1
run1d: meshfile = diode1D.msh
run1d: format := gnuplot
run1d: diode1D.msh

run2d: DIM := 2
run2d: meshfile = diode2D_mdr.grd
run2d: format := gmv

$(targets): run

run:
	@echo "Running example in $(DIM)D..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@perl -i -p -e 's/^\s*dimension\s*=.*/  dimension = $(DIM)/' $(infile)
	@perl -i -p -e 's/^\s*meshfile\s*=.*/  meshfile = $(meshfile)/' $(infile)
	@perl -i -p -e 's/^\s*resultpath\s*=.*/  resultpath = output_$(DIM)D/' $(infile)
	@perl -i -p -e 's/^\s*output_format\s*=.*/  output_format = $(format)/' $(infile)
	@./tibercad $(infile)


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output_*D tibercad


include $(topdir)/Make.rules
