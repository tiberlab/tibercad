# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

# default dimension of the example
DIM := 1

# default meshfile
meshfile := diode.msh
format := msh

# order of the mesh elements (1 or 2)
ORDER := 1


targets := run1d run2d

.PHONY: all clean run


all: main

main: $(objects)
	@echo "Linking main..."
	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) -o main $(LIBS) $(LDFLAGS) $(EXTERNAL_FLAGS)

run1d: DIM := 1
run1d: meshfile := diode.msh
run1d: format := msh
run2d: DIM := 2
run2d: meshfile := diode_mdr.grd
run2d: format := ise
$(targets): $(meshfile) run

run: all
	@echo "Running example in $(DIM)D..."
	@[ -d output ] || mkdir output
	@perl -i -p -e 's/^\s*dimension\s*=.*/dimension = $(DIM)/' diode.in
	@perl -i -p -e 's/^\s*meshfile\s*=.*/meshfile = $(meshfile)/' diode.in
	@perl -i -p -e 's/^\s*format\s*=.*/format = $(format)/' diode.in
	./main


clean:
	@rm -f *.o *~ main *.xda *.xta

distclean: clean
	@rm -f *.msh *.gmv* 
	@rm -rf output


include $(topdir)/Make.rules
