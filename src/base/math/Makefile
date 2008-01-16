# $Id$

topdir := ../../..

include $(topdir)/Make.common


#
# source files
fsrcfiles 	:= $(wildcard *.f)
csrcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.f, %.$(obj-suffix), $(fsrcfiles))
objects		+= $(patsubst %.C, %.$(obj-suffix), $(fsrcfiles))


all: $(objects)


include $(topdir)/Make.rules
