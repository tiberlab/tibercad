# $Id$

topdir := ../../..

include $(topdir)/Make.common


#
# source files
fsrcfiles := $(wildcard *.f)
csrcfiles := $(wildcard *.C)

#
# object files
target	:= $(patsubst %.f, %$(obj-suffix), $(fsrcfiles))
target	+= $(patsubst %.C, %$(obj-suffix), $(csrcfiles))



include $(topdir)/Make.rules

