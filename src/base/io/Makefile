# $Id$

topdir := ../../..

subdirs := Read_ISE_grid 

include $(topdir)/Make.common

csrcfiles := $(wildcard *.c)
cobjects := $(patsubst %.c, .$(ARCH)/%$(obj-suffix), $(csrcfiles))
target += $(cobjects)

include $(topdir)/Make.rules
