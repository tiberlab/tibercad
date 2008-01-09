# $Id$
#

topdir := ../../..

include $(topdir)/Make.common

# we have to compile each module
subdirs := petsc
ifeq ($(enable_pardiso), yes)
 subdirs += pardiso
endif


# source files
srcfiles 	:= $(wildcard *.C)

# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

.PHONY: all
all: $(objects)

echo:
	@echo $(subdirs)


include $(topdir)/Make.rules

