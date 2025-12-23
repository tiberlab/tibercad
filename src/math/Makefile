# $Id$

topdir := ../..

include $(topdir)/Make.common


#
# source files
fsrcfiles := $(wildcard *.f)
csrcfiles := $(wildcard *.cpp)

#
# object files
target	:= $(patsubst %.f, .$(ARCH)/%$(obj-suffix), $(fsrcfiles))
target	+= $(patsubst %.cpp, .$(ARCH)/%$(obj-suffix), $(csrcfiles))



include $(topdir)/Make.rules

