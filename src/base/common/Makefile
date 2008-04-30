# $Id$

topdir := ../../..

include $(topdir)/Make.common


#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))



all: $(objects)


TiberCad.$(obj-suffix): svnrevision.h

.PHONY: svnrevision.h
svnrevision.h:
	@echo "#define SVNREVISION $(SVNREVISION)" > svnrevision.h



include $(topdir)/Make.rules
