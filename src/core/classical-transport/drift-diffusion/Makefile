# $Id$

topdir := ../../../..

# subdirectories
subdirs		:= models


include $(topdir)/Make.common

# source files
srcfiles 	:= $(wildcard *.C)

# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

ifeq ($(build_modules), yes)
 target := shared
else
 target := $(objects)
endif


.PHONY: all
all: $(target)


shared: $(objects)
	@soname=lib`grep TIBER_MODULE DriftDiffusion.C | \
		sed 's/[[:space:]]*TIBER_MODULE([[:space:]]*[a-zA-Z]\+[[:space:]]*,[[:space:]]*\([a-zA-Z_]\+\)[[:space:]]*)/\1/'`.so ; \
		$(CXX) -shared $(CXXFLAGS) -Wl,-soname,$$soname \
			-o $(topdir)/lib/tibermodels/$$soname $(objects)


	

include $(topdir)/Make.rules
