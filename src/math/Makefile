# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
fsrcfiles := $(wildcard *.f)
csrcfiles := $(wildcard *.C)

#
# object files
target	:= $(patsubst %.f, .$(ARCH)/%$(obj-suffix), $(fsrcfiles))
target	+= $(patsubst %.C, .$(ARCH)/%$(obj-suffix), $(csrcfiles))



include $(topdir)/Make.rules

