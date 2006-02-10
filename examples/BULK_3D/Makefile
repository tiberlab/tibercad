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
DIM := 3

LIBS += $(topdir)/lib/libtensor.a


.PHONY: all clean run


all: main

main: $(objects)
	@echo "Linking main..."
	@$(CXX) $(CXXFLAGS) $(objects) $(TIBER_OBJECTS) -o main $(LIBS) $(LDFLAGS) 

run: all bulk.msh
	@echo "Running example..."
	./main


clean:
	@rm -f *.o *~ main *.xda *.xta *.msh


include $(topdir)/Make.rules
