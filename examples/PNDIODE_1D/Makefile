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
DIM := 1
# order of the mesh elements (1 or 2)
ORDER := 1


.PHONY: all clean run


all: main

main: $(objects)
	@echo "Linking main..."
	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) -o main $(LIBS) $(LDFLAGS) $(EXTERNAL_FLAGS)


run: all diode.msh
	@echo "Running example..."
	@[ -d output ] || mkdir output
	./main


clean:
	@rm -f *.o *~ main *.xda *.xta

distclean: clean
	@rm -f *.msh *.gmv* 
	@rm -rf output


include $(topdir)/Make.rules
