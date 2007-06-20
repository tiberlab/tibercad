# $Id$
#

topdir := ../../..

# we have to compile each module
subdirs := petsc pardiso

include $(topdir)/Make.common

# source files
srcfiles 	:= $(wildcard *.C)

# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

.PHONY: all
all: $(objects)



include $(topdir)/Make.rules

