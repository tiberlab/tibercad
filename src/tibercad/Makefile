# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))


all: tibercad

tibercad: $(objects)
	@echo "Compiling tibercad executable ..."
	@$(CXX) $(CXXFLAGS) $(objects) -o $(topdir)/bin/tibercad \
		-Wl,-export-dynamic -L$(topdir)/lib \
		-Wl,-rpath,$(topdir)/lib -ltibercad



include $(topdir)/Make.rules
