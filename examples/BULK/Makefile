# $Id$

topdir := ../..

include $(topdir)/Make.common

targets := run1d run2d run3d

.PHONY: run

DIM := 2

all: run

run1d: DIM := 1
run2d: DIM := 2
run3d: DIM := 3
$(targets): run%d: bulk%D.msh run
run2d_tri: run%d_tri: bulk2Dtri.msh 
	@echo "Running example in $(DIM)D..."
	@echo ""
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/  dimension = $(DIM)/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/  meshfile = bulk$(DIM)Dtri.msh/' bulk.in

run:
	@echo "Running example in $(DIM)D..."
	@echo ""
	@ln -f -s $(topdir)/bin/tibercad
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/  dimension = $(DIM)/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/  meshfile = bulk$(DIM)D.msh/' bulk.in
	@./tibercad bulk.in


clean:
	@rm -f *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output tibercad


include $(topdir)/Make.rules
