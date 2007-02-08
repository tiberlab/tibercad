# $Id$

topdir := ../..

include $(topdir)/Make.common


# default dimension of the example
DIM := 2

# default meshfile
meshfile := diode.msh
format := msh

# order of the mesh elements (1 or 2)
ORDER := 1


targets := run1d run2d

.PHONY: all clean run


all: run1d run2d

run1d: DIM := 1
run1d: meshfile := diode.msh
run2d: DIM := 2
run2d: meshfile := diode_mdr.grd
$(targets): $(meshfile) run

run: all
	@echo "Running example in $(DIM)D..."
	@ln -f -s $(topdir)/bin/tibercad
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/  dimension = $(DIM)/' diode.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/  meshfile = $(meshfile)/' diode.in
	@./tibercad diode.in


clean:
	@rm -f *~ *.xda *.xta

distclean: clean
	@rm -f *.msh *.gmv* 
	@rm -rf output tibercad


include $(topdir)/Make.rules
