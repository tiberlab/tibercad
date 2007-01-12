# $Id$

topdir := ../../..

subdirs := Read_ISE_grid

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))




all: $(objects)



include $(topdir)/Make.rules
