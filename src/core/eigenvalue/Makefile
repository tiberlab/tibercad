# $Id: Makefile 131 2006-11-29 12:04:01Z maufder $

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
