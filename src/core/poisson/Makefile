# $Id: Makefile 163 2007-01-12 12:37:28Z maufder $

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



all: $(objects)


include $(topdir)/Make.rules
