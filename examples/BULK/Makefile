# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

# dimension of the example
DIM := 2
# order of the mesh elements (1 or 2)
ORDER := 1

targets := run1d run2d run3d

.PHONY: run


all: main

main: $(objects)
	@echo "Linking main..."
	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) -o main $(LIBS) $(LDFLAGS) 

run1d: DIM := 1
run2d: DIM := 2
run3d: DIM := 3
$(targets): run%d: bulk%D.msh run
run2d_tri: run%d_tri: bulk2Dtri.msh main
	@echo "Running example in $(DIM)D..."
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/dimension = $(DIM)/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/meshfile = bulk$(DIM)Dtri.msh/' bulk.in
	./main

run: all
	@echo "Running example in $(DIM)D..."
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/dimension = $(DIM)/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/meshfile = bulk$(DIM)D.msh/' bulk.in
	./main


clean:
	@rm -f main *.xda *.xta *.ucd

distclean: clean
	@rm -f *.msh *.dat
	@rm -rf output


include $(topdir)/Make.rules
