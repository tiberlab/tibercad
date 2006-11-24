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

targets := run2d run3d

.PHONY: all clean run


all: main

main: $(objects)
	@echo "Linking main..."
	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) -o main $(LIBS) $(LDFLAGS) 

run2d: DIM := 2
run3d: DIM := 3
$(targets): run%d: bulk%D.msh run

run: all
	@echo "Running example in $(DIM)D..."
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/dimension = $(DIM)/' bulk.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/meshfile = bulk$(DIM)D.msh/' bulk.in
	./main


clean:
	@rm -f *.o *~ main *.xda *.xta

distclean: clean
	@rm -f *.msh
	@rm -rf output


include $(topdir)/Make.rules
