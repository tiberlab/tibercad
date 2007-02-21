# $Id$

topdir := ../..

include $(topdir)/Make.common

targets := run1d run2d run3d

.PHONY: run

DIM := 2

all: run

run1d: DIM := 1
run2d: DIM := 2
run2dtri: DIM := 2
run3d: DIM := 3
$(targets): run%d: bulk%D.msh run
run2d_tri: run%d_tri: bulk2Dtri.msh 
	@echo "Running example in 2D..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@perl -i -p -e 's/^\s*dimension\s*=.*/  dimension = 2/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/  meshfile = bulk2Dtri.msh/' bulk.in
	@perl -i -p -e 's/^\s*resultpath\s*=.*/  resultpath = output_tri/' bulk.in
	@./tibercad bulk.in

run:
	@echo "Running example in $(DIM)D..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@perl -i -p -e 's/^\s*dimension\s*=.*/  dimension = $(DIM)/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/  meshfile = bulk$(DIM)D.msh/' bulk.in
	@perl -i -p -e 's/^\s*resultpath\s*=.*/  resultpath = output_$(DIM)D/' bulk.in
	@./tibercad bulk.in


clean:
	@rm -f *~ *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output_*D output_tri tibercad


include $(topdir)/Make.rules
