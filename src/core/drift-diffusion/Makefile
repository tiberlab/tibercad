# $Id$

topdir := ../../..

include $(topdir)/Make.common

CXXFLAGS += -I/usr/pack/gsl-1.6-ma/include/


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



.PHONY: all clean


all: $(objects)

#-include $(srcfiles:.C=.d)


clean:
	@rm -f *.o *~ *.d


test:
	@$(CXX) $(CXXFLAGS) $(INCLUDE) -S DriftDiffusion.C
